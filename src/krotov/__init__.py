__version__ = '0.0.1'
import time
import numpy as np
import attr
import logging

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


@attr.s
class PulseOptions():
    """Options for the optimization of a control pulse

    Attributes:
        lambda_a (float): Krotov step size. This governs the overall magnitude
            of the pulse update. Large values result in small updates. Small
            values may lead to sharp spikes and numerical instability.
        shape (callable): Function S(t) in the range [0, 1] that scales the
            pulse update for the pulse value at t. This can be used to ensure
            boundary contitions (S(0) = S(T) = 0), and enforce smooth switch-on
            and switch-off
        filter (callable or None): A function that manipulates the pulse after
            each OCT iteration, e.g. by applying a spectral filter.
    """
    lambda_a = attr.ib()
    shape = attr.ib(default=lambda t: 1)
    filter = attr.ib(default=None)


class Result():
    """Result object for a Krotov optimization

    Attributes:
        objectives (list): A copy of the control objectives. Each item is an
            instance of :class:`Objective`, and we obtain a single
            set of controls that optimizes the average of all objectives.
        guess_controls (list): List of the original guess pulses
        optimized_controls (list): List of the optimized pulses, in the order
            corresponding to `guess_controls`
        tlist (numpy array): A copy of the time grid values
        tlist_midpoints (numpy array): points centered between the time points
            in `tgrid`
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
        self.objective = None
        self.tlist = tlist.copy()
        self.iters = objectives.copy()
        self.iter_seconds = []
        self.info_vals = []
        self.tau_vals = []
        self.guess_controls = guess_controls  # do not use copy
        self.optimized_controls = []
        self.updates = []
        self.states = [obj.target_state for obj in objectives]
        self.start_local_time = time.localtime()
        self.end_local_time = time.localtime()
        tlist_midpoints = []
        for i in range(len(tlist) - 1):
            tlist_midpoints.append(0.5 * (tlist[i+1] + tlist[i]))
        self.tlist_midpoints = np.array(tlist_midpoints)


def _extract_controls(objectives, pulse_options):
    """Extract a list of controls from the objectives (in arbitrary order), and
    check that all controls have proper options"""
    logger = logging.getLogger(__name__)
    controls = set()
    for i_obj, objective in enumerate(objectives):
        for i_ham, ham in enumerate(objective.H):
            if isinstance(ham, list):
                control = ham[1]
                if control in pulse_options:
                    controls.add(control)
                else:
                    raise ValueError(
                        "The control %s in the component %d of the "
                        "Hamiltonian of the objective %d does not have any "
                        "associated pulse options" % (control, i_ham, i_obj))
    if len(controls) != len(pulse_options):
        logger.warning(
            "pulse_options contains options for controls that are not in the "
            "objectives")
    return list(controls)


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
            There must be a mapping for each control.
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
    guess_controls = _extract_controls(objectives, pulse_options)
    result = Result(objectives, guess_controls, tlist)
    return result
