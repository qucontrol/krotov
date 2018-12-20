import logging
import time
import copy
from functools import partial
from typing import Optional

import numpy as np

from qutip import Qobj
from qutip.parallel import serial_map

from .result import Result
from .structural_conversions import (
    extract_controls, extract_controls_mapping, control_onto_interval,
    pulse_options_dict_to_list, pulse_onto_tlist,
    plug_in_pulse_values, discretize)
from .mu import derivative_wrt_pulse
from .info_hooks import chain

__all__ = ['optimize_pulses']


def _overlap(a, b) -> Optional[complex]:
    """Complex overlap of two quantum objects.

    If `a`, `b` are not quantum objects or are not compatible, return None.
    """
    if isinstance(a, Qobj) and isinstance(b, Qobj):
        if a.dims != b.dims:
            return None
        if a.type == b.type == 'oper':
            return complex((a.dag() * b).tr())
        else:
            return a.overlap(b)
    else:
        return None


def optimize_pulses(
    objectives,
    pulse_options,
    tlist,
    *,
    propagator,
    chi_constructor,
    mu=None,
    sigma=None,
    iter_start=0,
    iter_stop=5000,
    check_convergence=None,
    state_dependent_constraint=None,
    info_hook=None,
    modify_params_after_iter=None,
    storage='array',
    parallel_map=None,
    store_all_pulses=False
):
    r"""Use Krotov's method to optimize towards the given `objectives`

    Optimize all time-dependent controls found in the Hamiltonians of the given
    `objectives`.

    Args:
        objectives (list[Objective]): List of objectives
        pulse_options (dict): Mapping of time-dependent controls found in the
            Hamiltonians of the objectives to :class:`.PulseOptions` instances.
            There must be a mapping for each control. As numpy arrays are
            unhashable and thus cannot be used as dict keys, the options for a
            ``control`` that is an array must set using
            ``pulse_options[id(control)] = ...``
        tlist (numpy.ndarray): Array of time grid values, cf.
            :func:`~qutip.mesolve.mesolve`
        propagator (callable): Function that propagates the state backward or
            forwards in time by a single time step, between two points in
            `tlist`
        chi_constructor (callable): Function that calculates the boundary
            condition for the backward propagation. This is where the
            final-time functional (indirectly) enters the optimization.
        mu (None or callable): Function that calculates the derivative
            $\frac{\partial H}{\partial\epsilon}$ for an equation of motion
            $\dot{\phi}(t) = -i H[\phi(t)]$ of an abstract operator $H$ and an
            abstract state $\phi$. If None, defaults to
            :func:`krotov.mu.derivative_wrt_pulse`, which covers the standard
            Schrödinger and master equations. See :mod:`krotov.mu` for a
            full explanation of the role of `mu` in the optimization, and the
            required function signature.
        sigma (None or callable): Function that calculates the second-order
            Krotov term. If None, the first-order Krotov method is used.
        iter_start (int): The formal iteration number at which to start the
            optimization
        iter_stop (int): The iteration number after which to end the
            optimization, whether or not convergence has been reached
        check_convergence (None or callable): Function that determines whether
            the optimization has converged. If None, the optimization will only
            end when `iter_stop` is reached. See :mod:`krotov.convergence` for
            details.
        state_dependent_constraint (None or callable): Function that evaluates
            a state-dependent constraint. If None, optimize without any
            state-dependent constraint. Currently not implemented.
        info_hook (None or callable): Function that is called after each
            iteration of the optimization, for the purpose of analysis. Any
            value returned by `info_hook` (e.g. an evaluated functional J_T)
            will be stored, for each iteration, in the `info_vals` attribute of
            the returned :class:`.Result`. The `info_hook` must have the same
            signature as :func:`krotov.info_hooks.print_debug_information`. It
            should not modify its arguments except for `shared_data` in any
            way.
        modify_params_after_iter (None or callable): Function that is called
            after each iteration, which may modify its arguments for certain
            advanced use cases, such as dynamically adjusting `lambda_vals`, or
            applying spectral filters to the `optimized_pulses`. It has the
            same interface as `info_hook` but should not return anything. The
            `modify_params_after_iter` function is called immediately before
            `info_hook`, and can transfer arbitrary data to any subsequent
            `info_hook` via the `shared_data` argument.
        storage (callable): Storage constructor for the storage of
            propagated states. Must accept an integer parameter `N` and return
            an empty array of length `N`. The default value 'array' is
            equivalent to ``functools.partial(numpy.empty, dtype=object)``.
        parallel_map (callable or None): Parallel function evaluator. The
            argument must have the same specification as
            :func:`qutip.parallel.serial_map`,
            which is used when None is passed. Alternatives are
            :func:`qutip.parallel.parallel_map` or
            :func:`qutip.ipynbtools.parallel_map`.
        store_all_pulses (bool): Whether or not to store the optimized pulses
            from *all* iterations in :class:`.Result`.

    Returns:
        Result: The result of the optimization.

    Raises:
        ValueError: If any controls are not real-valued, or if any update
            shape is not a real-valued function in the range [0, 1].
    """
    logger = logging.getLogger('krotov')

    # Initialization
    logger.info("Initializing optimization with Krotov's method")
    if mu is None:
        mu = derivative_wrt_pulse
    if sigma is not None:
        raise NotImplementedError("Second order")
    if modify_params_after_iter is not None:
        # From a technical perspective, the `modify_params_after_iter` is
        # really just another info_hook, the only difference is the
        # convention that info_hooks shouldn't modify the parameters.
        if info_hook is None:
            info_hook = modify_params_after_iter
        else:
            info_hook = chain(modify_params_after_iter, info_hook)
    if state_dependent_constraint is not None:
        raise NotImplementedError("state_dependent_constraint")

    adjoint_objectives = [obj.adjoint for obj in objectives]
    if storage == 'array':
        storage = partial(np.empty, dtype=object)
    if parallel_map is None:
        parallel_map = serial_map

    (
        guess_controls,
        guess_pulses,
        pulses_mapping,
        lambda_vals,
        shape_arrays,
    ) = _initialize_krotov_controls(objectives, pulse_options, tlist)

    result = Result()
    result.start_local_time = time.localtime()

    # Initial forward-propagation
    tic = time.time()
    forward_states = parallel_map(
        _forward_propagation,
        list(range(len(objectives))),
        (objectives, guess_pulses, pulses_mapping, tlist, propagator, storage),
    )
    toc = time.time()

    fw_states_T = [states[-1] for states in forward_states]
    tau_vals = [
        _overlap(state_T, obj.target)
        for (state_T, obj) in zip(fw_states_T, objectives)
    ]

    info = info_hook(
        objectives=objectives,
        adjoint_objectives=adjoint_objectives,
        backward_states=None,
        forward_states=forward_states,
        optimized_pulses=guess_pulses,
        lambda_vals=lambda_vals,
        shape_arrays=shape_arrays,
        fw_states_T=fw_states_T,
        tau_vals=tau_vals,
        start_time=tic,
        stop_time=toc,
        iteration=0,
        shared_data={},
    )

    # Initialize result
    result.tlist = tlist
    result.objectives = objectives
    result.guess_controls = guess_controls
    result.controls_mapping = pulses_mapping
    result.info_vals.append(info)
    result.iters.append(0)
    result.tau_vals.append(tau_vals)
    if store_all_pulses:
        result.all_pulses.append(guess_pulses)
    result.states = fw_states_T

    for krotov_iteration in range(iter_start + 1, iter_stop + 1):

        logger.info("Started Krotov iteration %d" % krotov_iteration)
        tic = time.time()

        # Boundary condition for the backward propagation
        # -- this is where the functional enters the optimizaton
        chi_states = chi_constructor(fw_states_T, objectives, result.tau_vals)
        chi_norms = [chi.norm() for chi in chi_states]
        # normalizing χ improves numerical stability; the norm then has to be
        # taken into account when calculating Δϵ
        chi_states = [chi / nrm for (chi, nrm) in zip(chi_states, chi_norms)]

        # Backward propagation
        backward_states = parallel_map(
            _backward_propagation,
            list(range(len(objectives))),
            (
                chi_states,
                adjoint_objectives,
                guess_pulses,
                pulses_mapping,
                tlist,
                propagator,
                storage,
            ),
        )

        # Forward propagation and pulse update
        logger.info("Started forward propagation/pulse update")
        # forward_states_from_previous_iter = forward_states
        forward_states = [storage(len(tlist)) for _ in range(len(objectives))]
        for i_obj in range(len(objectives)):
            forward_states[i_obj][0] = objectives[i_obj].initial_state
        delta_eps = [
            np.zeros(len(tlist) - 1, dtype=np.complex128) for _ in guess_pulses
        ]
        optimized_pulses = copy.deepcopy(guess_pulses)
        for time_index in range(len(tlist) - 1):  # iterate over time intervals
            # pulse update
            for (i_pulse, guess_pulse) in enumerate(guess_pulses):
                for (i_obj, objective) in enumerate(objectives):
                    χ = backward_states[i_obj][time_index]
                    μ = mu(
                        objectives,
                        i_obj,
                        guess_pulses,
                        pulses_mapping,
                        i_pulse,
                        time_index,
                    )
                    Ψ = forward_states[i_obj][time_index]
                    update = _overlap(χ, μ(Ψ))  # ⟨χ|μ|Ψ⟩
                    update *= chi_norms[i_obj]
                    delta_eps[i_pulse][time_index] += update
                λₐ = lambda_vals[i_pulse]
                S_t = shape_arrays[i_pulse][time_index]
                Δϵ = (S_t / λₐ) * delta_eps[i_pulse][time_index].imag
                optimized_pulses[i_pulse][time_index] += Δϵ
            # forward propagation
            fw_states = parallel_map(
                _forward_propagation_step,
                list(range(len(objectives))),
                (
                    forward_states,
                    objectives,
                    optimized_pulses,
                    pulses_mapping,
                    tlist,
                    time_index,
                    propagator,
                ),
            )
            # storage
            for i_obj in range(len(objectives)):
                forward_states[i_obj][time_index + 1] = fw_states[i_obj]
        logger.info("Finished forward propagation/pulse update")
        fw_states_T = fw_states
        tau_vals = [
            _overlap(fw_state_T, obj.target)
            for (fw_state_T, obj) in zip(fw_states_T, objectives)
        ]

        toc = time.time()

        # Display information about iteration
        if info_hook is None:
            info = None
        else:
            info = info_hook(
                objectives=objectives,
                adjoint_objectives=adjoint_objectives,
                backward_states=backward_states,
                forward_states=forward_states,
                fw_states_T=fw_states_T,
                optimized_pulses=optimized_pulses,
                lambda_vals=lambda_vals,
                shape_arrays=shape_arrays,
                tau_vals=tau_vals,
                start_time=tic,
                stop_time=toc,
                shared_data={},
                iteration=krotov_iteration,
            )
        # Update optimization `result` with info from finished iteration
        result.iters.append(krotov_iteration)
        result.iter_seconds.append(int(toc - tic))
        result.info_vals.append(info)
        result.tau_vals.append(tau_vals)
        result.optimized_controls = optimized_pulses
        if store_all_pulses:
            result.all_pulses.append(optimized_pulses)
        result.states = fw_states_T

        logger.info("Finished Krotov iteration %d" % krotov_iteration)

        # Convergence check
        msg = None
        if check_convergence is not None:
            msg = check_convergence(result)
        if bool(msg) is True:  # this is not an anti-pattern!
            result.message = "Reached convergence"
            if isinstance(msg, str):
                result.message += ": " + msg
            break
        else:
            # prepare for next iteration
            guess_pulses = optimized_pulses

    else:  # optimization finished without `check_convergence` break

        result.message = "Reached %d iterations" % iter_stop

    # Finalize
    result.end_local_time = time.localtime()
    for i, pulse in enumerate(optimized_pulses):
        result.optimized_controls[i] = pulse_onto_tlist(pulse)
    return result


def _initialize_krotov_controls(objectives, pulse_options, tlist):
    """Extract discretized guess controls and pulses from `objectives`, and
    return them with the associated mapping and option data"""
    guess_controls = extract_controls(objectives)
    pulses_mapping = extract_controls_mapping(objectives, guess_controls)
    options_list = pulse_options_dict_to_list(pulse_options, guess_controls)
    guess_controls = [discretize(control, tlist) for control in guess_controls]
    for control in guess_controls:
        if np.iscomplexobj(control):
            raise ValueError(
                "All controls must be real-valued. Complex controls must be "
                "split into an independent real and imaginary part in the "
                "objectives before passing them to the optimization")
    guess_pulses = [  # defined on the tlist intervals
        control_onto_interval(control)
        for control in guess_controls]
    lambda_vals = [options.lambda_a for options in options_list]
    shape_arrays = []
    for options in options_list:
        S = discretize(options.shape, tlist, args=())
        shape_arrays.append(control_onto_interval(S))
    for shape_array in shape_arrays:
        if np.iscomplexobj(shape_array):
            raise ValueError(
                "Update shapes (shape attribute in PulseOptions) must be "
                "real-valued")
        if np.min(shape_array) < 0 or np.max(shape_array) > 1.01:
            # 1.01 accounts for rounding errors: In principle, shapes > 1 are
            # not a problem, but then it cancels with λₐ, which makes things
            # unnecessarily confusing.
            raise ValueError(
                "Update shapes (shape attribute in PulseOptions) must have "
                "values in the range [0, 1]")
    return (
        guess_controls, guess_pulses, pulses_mapping, lambda_vals,
        shape_arrays)


def _forward_propagation(
    i_objective,
    objectives,
    pulses,
    pulses_mapping,
    tlist,
    propagator,
    storage,
    store_all=True,
):
    """Forward propagation of the initial state of a single objective over the
    entire `tlist`"""
    logger = logging.getLogger('krotov')
    logger.info(
        "Started initial forward propagation of objective %d", i_objective
    )
    obj = objectives[i_objective]
    state = obj.initial_state
    if store_all:
        storage_array = storage(len(tlist))
        storage_array[0] = state
    mapping = pulses_mapping[i_objective]
    for time_index in range(len(tlist) - 1):  # index over intervals
        H = plug_in_pulse_values(obj.H, pulses, mapping[0], time_index)
        c_ops = [
            plug_in_pulse_values(c_op, pulses, mapping[ic + 1], time_index)
            for (ic, c_op) in enumerate(obj.c_ops)
        ]
        dt = tlist[time_index + 1] - tlist[time_index]
        state = propagator(H, state, dt, c_ops)
        if store_all:
            storage_array[time_index + 1] = state
    logger.info(
        "Finished initial forward propagation of objective %d", i_objective
    )
    if store_all:
        return storage_array
    else:
        return state


def _backward_propagation(
    i_state,
    chi_states,
    adjoint_objectives,
    pulses,
    pulses_mapping,
    tlist,
    propagator,
    storage,
):
    """Backward propagation of chi_states[i_state] over the entire `tlist`"""
    logger = logging.getLogger('krotov')
    logger.info("Started backward propagation of state %d", i_state)
    state = chi_states[i_state]
    obj = adjoint_objectives[i_state]
    storage_array = storage(len(tlist))
    storage_array[-1] = state
    mapping = pulses_mapping[i_state]
    for time_index in range(len(tlist) - 2, -1, -1):  # index bw over intervals
        H = plug_in_pulse_values(
            obj.H, pulses, mapping[0], time_index, conjugate=True
        )
        c_ops = [
            plug_in_pulse_values(c_op, pulses, mapping[ic + 1], time_index)
            for (ic, c_op) in enumerate(obj.c_ops)
        ]
        dt = tlist[time_index + 1] - tlist[time_index]
        state = propagator(H, state, dt, c_ops, backwards=True)
        storage_array[time_index] = state
    logger.info("Finished backward propagation of state %d", i_state)
    return storage_array


def _forward_propagation_step(
    i_state,
    states,
    objectives,
    pulses,
    pulses_mapping,
    tlist,
    time_index,
    propagator,
):
    """Forward-propagate states[i_state] by a single time step"""
    state = states[i_state][time_index]
    obj = objectives[i_state]
    mapping = pulses_mapping[i_state]
    H = plug_in_pulse_values(obj.H, pulses, mapping[0], time_index)
    c_ops = [
        plug_in_pulse_values(c_op, pulses, mapping[ic + 1], time_index)
        for (ic, c_op) in enumerate(obj.c_ops)
    ]
    dt = tlist[time_index + 1] - tlist[time_index]
    return propagator(H, state, dt, c_ops)
