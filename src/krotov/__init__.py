__version__ = '0.0.1'
import time
import numpy as np
import attr


__all__ = ['Result', 'Objective', 'optimize_pulse']


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


class Result():
    """Result object for a Krotov optimization

    Attributes:
        objectives (list): A copy of the control objectives. Each item is an
            instance of :class:`KrotovObjective`, and we obtain a single
            set of controls that optimizes the average of all objectives.
        tlist (numpy array): A copy of the time grid values
        control_tlist (numpy array): The time grid for the `controls`, i.e.
            the points centered between the time points in `tgrid`
        iters (list of int): Iteration numbers
        iter_seconds (list of int): for each iteration number, the number of
            seconds that were spent in the optimization
        info_vals (list): For each iteration, the return value of `info_hook`,
            or None
        tau_vals (list of list): for each iteration, a list of complex overlaps
            between the forward-propagated states and the target states for
            each objective.
        controls (list of list): If the propagation was performed with
            ``store_all=True``, for each iteration, a list of optimized control
            fields (as an array). If ``store_all=False``, a list containing a
            single element, the list of final optimized control fields.
        states (list): for each objective, a list of states
            (:class:`qutip.Qobj` instances) for each value in
            `tlist`, obtained from propagation under the final optimized
            control fields.
        start_local_time (time.struct_time): Time stamp of when the
            optimization started
        end_local_time (time.struct_time): Time stamp of when the optimization
            ended
    """

    def __init__(self, objectives, controls, tlist):
        self.objective = None
        self.tlist = tlist.copy()
        self.iters = objectives.copy()
        self.iter_seconds = []
        self.info_vals = []
        self.tau_vals = []
        self.controls = [controls.copy()]
        self.states = [obj.target_state for obj in objectives]
        self.start_local_time = time.localtime()
        self.end_local_time = time.localtime()
        control_tlist = []
        for i in range(len(tlist) - 1):
            control_tlist.append(0.5 * (tlist[i+1] + tlist[i]))
        self.control_tlist = np.array(control_tlist)


def optimize_pulse(
        objectives, controls, tlist, propagator, chi_constructor,
        sigma=None, iter_start=0, iter_stop=5000, check_convergence=None,
        state_dependent_constraint=None, info_hook=None, store_all=False):
    """Use Krotov's method to optimize `controls` towards the given
    `objectives`

    Args:
        objectives (list): List of objectives
        controls (list): List of time-dependent functions that appear in the
            Hamiltonians of the `objectives` (see time-dependent operators for
            :func:`~qutip.mesolve.mesolve`)
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
            :class:`KrotovResult`.
        store_all (bool): Whether or not to store the optimized control fields
            from *all* iterations in :class:`KrotovResult`. If False,
            :class:`KrotovResult` will only include information from the final
            iteration

    Returns:
        KrotovResult: The result of the optimization.
    """
    result = Result(objectives, controls, tlist)
    return result
