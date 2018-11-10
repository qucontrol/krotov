import copy
import time

from .structural_conversions import _tlist_midpoints

__all__ = ['Result']


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
