import time
from textwrap import dedent

from .objective import Objective
from .structural_conversions import (
    _nested_list_shallow_copy)

__all__ = ['Result']


class Result():
    """Result object for a Krotov optimization

    Attributes:
        objectives (list): The control objectives
        tlist (numpy array): The time grid values
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
    time_fmt = "%Y-%m-%d %H:%M:%S"

    def __init__(self):
        self.objectives = []
        self.tlist = []
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

    def __str__(self):
        return dedent(r'''
        Krotov Optimization Result
        --------------------------
        - Started at {start_local_time}
        - Number of objectives: {n_objectives}
        - Number of iterations: {n_iters}
        - Ended at {end_local_time}
        '''.format(
            start_local_time=self.start_local_time_str,
            n_objectives=len(self.objectives),
            n_iters=len(self.iters)-1,  # do not count zero iteration
            end_local_time=self.end_local_time_str,
        )).strip()

    def __repr__(self):
        return self.__str__()

    @property
    def start_local_time_str(self):
        """The `start_local_time` attribute formatted as a string"""
        if self.start_local_time is not None:
            return time.strftime(self.time_fmt, self.start_local_time)
        else:
            return 'n/a'

    @property
    def end_local_time_str(self):
        """The `end_local_time` attribute formatted as a string"""
        if self.end_local_time is not None:
            return time.strftime(self.time_fmt, self.end_local_time)
        else:
            return 'n/a'

    @property
    def optimized_objectives(self):
        """A copy of the objectives with the `optimized_controls` plugged in"""
        objectives = []
        for (i, obj) in enumerate(self.objectives):
            H = _plug_in_optimized_controls(
                obj.H, self.optimized_controls, self.controls_mapping[i][0])
            c_ops = [
                _plug_in_optimized_controls(
                    c_op, self.optimized_controls,
                    self.controls_mapping[i][j+1])
                for (j, c_op) in enumerate(obj.c_ops)]
            objectives.append(
                Objective(
                    H=H,
                    initial_state=obj.initial_state,
                    target_state=obj.target_state,
                    c_ops=c_ops))
        return objectives


def _plug_in_optimized_controls(H, controls, mapping):
    """Auxilliary routine to :attr:`Result.optimized_objectives`"""
    H = _nested_list_shallow_copy(H)
    for (control, control_mapping) in zip(controls, mapping):
        for i in control_mapping:
            H[i][1] = control
    return H
