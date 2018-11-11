import copy
import time

__all__ = ['Result']


class Result():
    """Result object for a Krotov optimization

    Attributes:
        objectives (list): The control objectives, with "extracted" controls.
        tlist (numpy array): The time grid values
        tlist_midpoints (numpy array): Points centered between the time points
            in `tlist`
        iters (list of int): Iteration numbers, starting at 0.
        iter_seconds (list of int): for each iteration number, the number of
            seconds that were spent in the optimization
        info_vals (list): For each iteration, the return value of `info_hook`,
            or None
        tau_vals (list of list): for each iteration, a list of complex overlaps
            between the forward-propagated state and the target state for
            each objective.
        guess_controls (list): List of the guess controls in array format
        optimized_controls (list): List of the optimized control fields, in the
            order corresponding to `guess_controls`
        controls_mapping (list): A nested list that indicates where in
            `objectives` the `guess_controls` and `optimized_controls` are used
            (as returned by :func:`.extract_controls_mapping`)
        all_pulses (list of list): If the optimization was performed with
            ``store_all_pulses=True``, for each iteration, a list of the
            optimized pulses (in the order corresponding to `guess_controls`).
            These pulses are defined at midpoints of the `tlist` intervals.
            Empty list if ``store_all_pulses=False``
        states (list): for each objective, a list of states
            (:class:`qutip.Qobj` instances) for each value in
            `tlist`, obtained from propagation under the final optimized
            control fields.
        start_local_time (time.struct_time): Time stamp of when the
            optimization started
        end_local_time (time.struct_time): Time stamp of when the optimization
            ended
    """

    def __init__(self):
        self.objectives = []
        self.tlist = []
        self.tlist_midpoints = []
        self.iters = []
        self.iter_seconds = []
        self.info_vals = []
        self.tau_vals = []
        self.guess_controls = []
        self.optimized_controls = []
        self.controls_mapping = []
        self.all_pulses = []
        self.states = []
        self.start_local_time = None
        self.end_local_time = None
