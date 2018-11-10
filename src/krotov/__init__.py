__version__ = '0.0.1'
import time
import numpy as np
import attr
import copy
import logging
from collections import OrderedDict, defaultdict

# expose submodules
from . import shapes
from . import pulse_conversion


__all__ = ['Result', 'Objective', 'PulseOptions', 'optimize_pulses']


@attr.s
class Objective():
    """A single objective for optimization with Krotov's method

    Attributes:
        initial_state (qutip.Qobj): The initial state
        target_state (qutip.Qobj): The desired target state
        H (list): The time-dependent Hamiltonian,
            cf. :func:`qutip.mesolve.mesolve`. This includes the control
            fields.
        c_ops (list):  List of collapse operators,
            cf. :func:`~qutip.mesolve.mesolve`.
    """
    H = attr.ib()
    initial_state = attr.ib()
    target_state = attr.ib()
    c_ops = attr.ib(default=[])

    def __copy__(self):
        # When we use copy.copy(objective), we want a
        # semi-deep copy where nested lists in the Hamiltonian and the c_ops
        # are re-created (copy by value), but non-list elements are copied by
        # reference. The _extract_controls relies heavily on this exact
        # behavior.
        return Objective(
            H=_nested_list_shallow_copy(self.H),
            initial_state=self.initial_state,
            target_state=self.target_state,
            c_ops=[_nested_list_shallow_copy(c) for c in self.c_ops])


def _nested_list_shallow_copy(l):
    if isinstance(l, list):
        return [copy.copy(h) if isinstance(h, list) else h for h in l]
    else:
        return l


@attr.s
class PulseOptions():
    """Options for the optimization of a control pulse

    Attributes:
        lambda_a (float): Krotov step size. This governs the overall magnitude
            of the pulse update. Large values result in small updates. Small
            values may lead to sharp spikes and numerical instability.
        shape (callable): Function S(t) in the range [0, 1] that scales the
            pulse update for the pulse value at t. This can be used to ensure
            boundary conditions (S(0) = S(T) = 0), and enforce smooth switch-on
            and switch-off
        filter (callable or None): A function that manipulates the pulse after
            each OCT iteration, e.g. by applying a spectral filter.
    """
    lambda_a = attr.ib()
    shape = attr.ib(default=lambda t: 1)
    filter = attr.ib(default=None)


def _tlist_midpoints(tlist):
    """Calculate array of midpoints in `tlist`"""
    tlist_midpoints = []
    for i in range(len(tlist) - 1):
        tlist_midpoints.append(0.5 * (tlist[i+1] + tlist[i]))
    return np.array(tlist_midpoints)


class Result():
    """Result object for a Krotov optimization

    Attributes:
        objectives (list): A copy of the control objectives. Each item is an
            instance of :class:`Objective`, and we obtain a single
            set of controls that optimizes the average of all objectives.
            The `objectives` will be in "extracted" form: time-dependent
            Hamiltonian components, which normally in QuTiP are in the form of
            a list ``[H, control]``, will instead be in the form ``[H, i]``,
            where ``i`` is the index of the corresponding control in
            `guess_controls` or `optimized_controls`.
        guess_controls (list): List of the original guess pulses
        optimized_controls (list): List of the optimized control fields, in the
            order corresponding to `guess_controls`
        tlist (numpy array): A copy of the time grid values
        tlist_midpoints (numpy array): points centered between the time points
            in `tlist`
        iters (list of int): Iteration numbers
        iter_seconds (list of int): for each iteration number, the number of
            seconds that were spent in the optimization
        info_vals (list): For each iteration, the return value of `info_hook`,
            or None
        tau_vals (list of list): for each iteration, a list of complex overlaps
            between the forward-propagated states and the target states for
            each objective.
        updates (list of list): If the optimization was performed with
            ``store_updates=True``, for each iteration, a list of the updates
            for all control fields (in the order corresponding to
            `guess_controls`). These updates are defined at midpoints of the
            `tlist` intervals. Empty list if ``store_updates=True``
        states (list): for each objective, a list of states
            (:class:`qutip.Qobj` instances) for each value in
            `tlist`, obtained from propagation under the final optimized
            control fields.
        start_local_time (time.struct_time): Time stamp of when the
            optimization started
        end_local_time (time.struct_time): Time stamp of when the optimization
            ended
    """

    def __init__(self, objectives, guess_controls, tlist):
        self.objectives = [copy.copy(o) for o in objectives]
        self.tlist = tlist.copy()
        self.iters = []
        self.iter_seconds = []
        self.info_vals = []
        self.tau_vals = []
        self.guess_controls = guess_controls  # do not use copy
        self.optimized_controls = []
        self.updates = []
        self.states = [obj.target_state for obj in objectives]
        self.start_local_time = time.localtime()
        self.end_local_time = time.localtime()
        self.tlist_midpoints = _tlist_midpoints(tlist)


def _find_in_list(val, list_to_search):
    """Return index of `val` in `list_to_search`, or -1

    Works even if `val` is a numpy array. In this case, comparison is by object
    identity.
    """
    if isinstance(val, np.ndarray):
        for i, v in enumerate(list_to_search):
            if v is val:
                return i
        return -1
    else:
        try:
            return list_to_search.index(val)
        except ValueError:
            return -1


def _extract_controls(objectives, pulse_options):
    """Extract a list of controls from the `objectives`, modify the
    Hamiltonians into an "extracted" form, and check that all controls have
    proper options

    Args:
        objectives (list): List of :class:`Objective` instances
        pulse_options (dict): Mapping of time-dependent controls found in the
            Hamiltonians of the objectives to :class:`PulseOptions` instances.
            There must be a mapping for each control.

    Returns:
        tuple: A tuple of the following:
        - list of `controls` extracted from `objective`
        - list of "pointers" for where each control occurs in the `objectives`.
          Each element is a list of tuples (objective-index,
          ham-component-index), where the objective-index gives the index of an
          objective that contains the control, and ham-component-index gives
          the index of a component of the Hamiltonian that is linear in the
          control
        - list of values from `pulse_options`, in the same order as `controls`
        - modified copy of `objectives`, where for each Hamiltonian component
          that contains a control, that control has been replaced with an
          integer that specifies in index of the corresponding control in
          `controls`.

    Example:

        >>> import qutip
        >>> X, Y, Z = qutip.Qobj(), qutip.Qobj(), qutip.Qobj() # dummy Hams
        >>> u1, u2 = np.array([]), np.array([])                # dummy controls
        >>> psi0, psi_tgt = qutip.Qobj(), qutip.Qobj()         # dummy states

        >>> H1 = [X, [Y, u1], [Z, u2]]  # ham for first objective
        >>> H2 = [X, [Y, u2]]           # ham for second objective
        >>> objectives = [
        ...     Objective(H1, psi0, psi_tgt),
        ...     Objective(H2, psi0, psi_tgt)]
        >>> pulse_options = {
        ...     id(u1): PulseOptions(lambda_a=1.0),
        ...     id(u2): PulseOptions(lambda_a=1.0)}

        >>> controls, control_map, options, objectives = _extract_controls(
        ...     objectives, pulse_options)
        >>> assert controls == [u1, u2]
        >>> assert objectives[0].H == [X, [Y, 0], [Z, 1]]
        >>> assert objectives[1].H == [X, [Y, 1]]
        >>> assert control_map[0] == [(0, 1)]          # where does u1 occur?
        >>> assert control_map[1] == [(0, 2), (1, 1)]  # where does u2 occur?
        >>> assert options[0] == pulse_options[id(u1)]
        >>> assert options[1] == pulse_options[id(u2)]
    """
    logger = logging.getLogger(__name__)
    controls = []
    controls_map = []
    options_list = []
    objectives = [copy.copy(o) for o in objectives]
    for i_obj, objective in enumerate(objectives):
        for i_ham, ham in enumerate(objective.H):
            if isinstance(ham, list):
                assert len(ham) == 2
                control = ham[1]
                i_control = _find_in_list(control, controls)
                if i_control >= 0:
                    controls_map[i_control].append((i_obj, i_ham))
                else:  # this is a control we haven't seen before
                    try:
                        try:
                            options_list.append(pulse_options[control])
                        except TypeError:  # control is numpy array
                            options_list.append(pulse_options[id(control)])
                    except KeyError:
                        raise ValueError(
                            "The control %s in the component %d of the "
                            "Hamiltonian of the objective %d does not have "
                            "any associated pulse options"
                            % (control, i_ham, i_obj))
                    controls.append(control)
                    controls_map.append([(i_obj, i_ham)])
                    i_control = len(controls) - 1
                ham[1] = i_control
    if len(controls) != len(pulse_options):
        logger.warning(
            "pulse_options contains options for controls that are not in the "
            "objectives")
    return controls, controls_map, options_list, objectives


def optimize_pulses(
        objectives, pulse_options, tlist, propagator, chi_constructor,
        sigma=None, iter_start=0, iter_stop=5000, check_convergence=None,
        state_dependent_constraint=None, info_hook=None, store_updates=False):
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
            forwards in time by a single time step, between to points in
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
        store_updates (bool): Whether or not to store the pulse updates
            from *all* iterations in :class:`Result`. These could be used to
            calculate the optimized pulses in each iteration.

    Returns:
        Result: The result of the optimization.
    """
    (guess_controls, control_mappings,
        options_list, objectives_in_extracted_form) = _extract_controls(
            objectives, pulse_options)
    guess_pulses = [  # defined on the tlist intervals
        pulse_conversion.control_onto_interval(
            control, tlist, _tlist_midpoints(tlist))
        for control in guess_controls]
    guess_controls = [  # convert guess controls to arrays, on tlist
        pulse_conversion.pulse_onto_tlist(pulse)
        for pulse in guess_pulses]
    result = Result(objectives_in_extracted_form, guess_controls, tlist)
    return result
