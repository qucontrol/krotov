import time

from .result import Result
from .structural_conversions import (
    extract_controls, control_onto_interval, pulse_onto_tlist,
    _tlist_midpoints)

__all__ = ['optimize_pulses']


def optimize_pulses(
        objectives, pulse_options, tlist, propagator, chi_constructor,
        sigma=None, iter_start=0, iter_stop=5000, check_convergence=None,
        state_dependent_constraint=None, info_hook=None, storage=list,
        store_all_pulses=False):
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
        storage (type or callable): Storage constructor for the storage of
            propagated states. This is the main memory requirement in the
            optimization, so this parameter gives you the option to provide an
            out-of-memory container.
        store_all_pulses (bool): Whether or not to store the optimized pulses
            from *all* iterations in :class:`Result`.

    Returns:
        Result: The result of the optimization.
    """
    (guess_controls, control_mappings,
        options_list, objectives_in_extracted_form) = extract_controls(
            objectives, pulse_options)
    tlist_midpoints = _tlist_midpoints(tlist)
    guess_pulses = [  # defined on the tlist intervals
        control_onto_interval(control, tlist, tlist_midpoints)
        for control in guess_controls]
    guess_controls = [  # convert guess controls to arrays, on tlist
        pulse_onto_tlist(pulse) for pulse in guess_pulses]

    result = Result()
    result.tlist = tlist
    result.tlist_midpoints = tlist_midpoints
    result.start_localtime = time.localtime()
    result.objectives = objectives_in_extracted_form
    result.guess_controls = guess_controls
    result.control_mappings = control_mappings
    if store_all_pulses:
        result.all_pulses.append(guess_pulses)

    # TODO: optimize

    result.end_local_time = time.localtime()
    return result
