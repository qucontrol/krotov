import logging
import time
from functools import partial

import numpy as np

from .result import Result
from .structural_conversions import (
    extract_controls, extract_controls_mapping, control_onto_interval,
    pulse_options_dict_to_list, pulse_onto_tlist, _tlist_midpoints,
    plug_in_pulse_values)

__all__ = ['optimize_pulses']


def serial_map(task, values, task_args=tuple(), task_kwargs=None, **kwargs):
    """Apply function `task` to all `values`

    Equivalent to::

        [task(value, *task_args, **task_kwargs) for value in values]
    """
    if task_kwargs is None:
        task_kwargs = {}
    return [task(value, *task_args, **task_kwargs) for value in values]


def _forward_propagation(
        i_objective, objectives, pulses, pulses_mapping, tlist,
        propagator, storage, store_all=True):
    """Do a forward propagation of the initial state of a single objective"""
    logger = logging.getLogger('krotov')
    logger.info(
        "Started initial forward propagation of objective %d", i_objective)
    obj = objectives[i_objective]
    state = obj.initial_state
    N = len(tlist) - 1  # number of time intervals
    if store_all:
        storage_array = storage(N)
    mapping = pulses_mapping[i_objective]
    for time_index in range(N):
        H = plug_in_pulse_values(obj.H, pulses, mapping[0], time_index)
        c_ops = [
            plug_in_pulse_values(c_op, pulses, mapping[ic+1], time_index)
            for (ic, c_op) in enumerate(obj.c_ops)]
        dt = tlist[time_index+1] - tlist[time_index]
        state = propagator(H, state, dt, c_ops)
        if store_all:
            storage_array[time_index] = state
    logger.info(
        "Finished initial forward propagation of objective %d", i_objective)
    if store_all:
        return storage_array
    else:
        return state


def optimize_pulses(
        objectives, pulse_options, tlist, propagator, chi_constructor,
        sigma=None, iter_start=0, iter_stop=5000, check_convergence=None,
        state_dependent_constraint=None, info_hook=None, storage='array',
        parallel_map=None, store_all_pulses=False):
    """Use Krotov's method to optimize towards the given `objectives`

    Optimize all time-dependent controls found in the Hamiltonians of the given
    `objectives`.

    Args:
        objectives (list): List of objectives
        pulse_options (dict): Mapping of time-dependent controls found in the
            Hamiltonians of the objectives to :class:`PulseOptions` instances.
            There must be a mapping for each control. As numpy arrays are
            unhashable and thus cannot be used as dict keys, the options for a
            ``control`` that is an array must set using
            ``pulse_options[id(control)] = ...``
        tlist (numpy array): Array of time grid values, cf.
            :func:`~qutip.mesolve.mesolve`
        propagator (callable): Function that propagates the state backward or
            forwards in time by a single time step, between two points in
            `tlist`
        chi_constructor (callable): Function that calculates the boundary
            condition for the backward propagation.
        sigma (None or callable): Function that calculates the second-order
            Krotov term. If None, the first-order Krotov method is used.
        iter_start (int): The formal iteration number at which to start the
            optimization
        iter_stop (int): The iteration number after which to end the
            optimization, whether or not convergence has been reached
        check_convergence (None or callable): Function that determines whether
            the optimization has converged. If None, the optimization will only
            end when `iter_stop` is reached.
        state_dependent_constraint (None or callable): Function that evaluates
            a state-dependent constraint. If None, optimize without any
            state-dependent constraint.
        info_hook (None or callable): Function that is called after each
            iteration of the optimization. Any value returned by `info_hook`
            (e.g. an evaluated functional J_T) will be stored, for each
            iteration, in the `info_vals` attribute of the returned
            :class:`Result`.
        storage (callable): Storage constructor for the storage of
            propagated states. Must accept an integer parameter `N` and return
            an empty array of length `N`. The default value 'array' is
            equivalent to ``functools.partial(numpy.empty, dtype=object)``.
        parallel_map (callable or None): Parallel function evaluator. The
            argument must have the same specification as :func:`serial_map`,
            which is used when None is passed. Alternatives are
            :func:`qutip.parallel.parallel_map` or
            :func:`qutip.ipynbtools.parallel_map`.
        store_all_pulses (bool): Whether or not to store the optimized pulses
            from *all* iterations in :class:`Result`.

    Returns:
        Result: The result of the optimization.
    """
    logger = logging.getLogger('krotov')

    # Initialization
    logger.info("Initializing optimization with Krotov's method")
    guess_controls = extract_controls(objectives)
    controls_mapping = extract_controls_mapping(objectives, guess_controls)
    options_list = pulse_options_dict_to_list(pulse_options, guess_controls)
    tlist_midpoints = _tlist_midpoints(tlist)
    guess_pulses = [  # defined on the tlist intervals
        control_onto_interval(control, tlist, tlist_midpoints)
        for control in guess_controls]
    guess_controls = [  # convert guess controls to arrays, on tlist
        pulse_onto_tlist(pulse) for pulse in guess_pulses]
    adjoint_objectives = [obj.adjoint for obj in objectives]
    if storage == 'array':
        storage = partial(np.empty, dtype=object)
    if parallel_map is None:
        parallel_map = serial_map

    M = len(objectives)  # number of objectives
    N = len(tlist) - 1   # number of time intervals

    result = Result()
    result.tlist = tlist
    result.tlist_midpoints = tlist_midpoints
    result.start_local_time = time.localtime()
    result.objectives = objectives
    result.guess_controls = guess_controls
    result.controls_mapping = controls_mapping

    # Initial forward-propagation
    tic = time.clock()
    forward_states = parallel_map(
        _forward_propagation, list(range(M)), (
            objectives, guess_controls, controls_mapping, tlist, propagator,
            storage))
    toc = time.clock()
    states_T = [states[-1] for states in forward_states]  # for each objective
    result.states = states_T
    result.iters.append(0)
    result.iter_seconds.append(int(toc - tic))
    result.tau_vals.append([
        states_T[i].overlap(objectives[i].target_state) for i in range(M)])
    if store_all_pulses:
        result.all_pulses.append(guess_pulses)

    # Boundary condition for the backward propagation
    # -- this is where the functional enters the optimizaton
    chi_states = chi_constructor(states_T, objectives, result.tau_vals)
    chi_norms = [chi.norm() for chi in chi_states]
    chi_states = [chi_states[i]/chi_norms[i] for i in range(M)]

    for krotov_iteration in range(iter_start, iter_stop+1):
        pass

    result.end_local_time = time.localtime()
    return result
