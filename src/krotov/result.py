import time
import pickle
import logging
from textwrap import dedent

from .objectives import Objective, _ControlPlaceholder
from .structural_conversions import _nested_list_shallow_copy

__all__ = ['Result']


class Result:
    """Result object for a Krotov optimization

    Attributes:
        objectives (list[Objective]): The control objectives
        tlist (numpy.ndarray): The time grid values
        iters (list[int]): Iteration numbers, starting at 0.
        iter_seconds (list[int]): for each iteration number, the number of
            seconds that were spent in the optimization
        info_vals (list): For each iteration, the return value of `info_hook`,
            or None
        tau_vals (list[list[complex]): for each iteration, a list of complex
            overlaps between the forward-propagated state and the target state
            for each objective, assuming :attr:`.Objective.target` contains the
            target state.
        guess_controls (list[numpy.ndarray]): List of the guess controls in
            array format
        optimized_controls (list[numpy.ndarray]): List of the optimized control
            fields, in the order corresponding to :attr:`guess_controls`
        controls_mapping (list): A nested list that indicates where in
            :attr:`objectives` the :attr:`guess_controls` and
            :attr:`optimized_controls` are used (as returned by
            :func:`.extract_controls_mapping`)
        all_pulses (list): If the optimization was performed with
            ``store_all_pulses=True``, for each iteration, a list of the
            optimized pulses (in the order corresponding to
            :attr:`guess_controls`).
            These pulses are defined at midpoints of the `tlist` intervals.
            Empty list if ``store_all_pulses=False``
        states (list[list[qutip.Qobj]]): for each objective, a list of states
            for each value in `tlist`, obtained from propagation under the
            final optimized control fields.
        start_local_time (time.struct_time): Time stamp of when the
            optimization started
        end_local_time (time.struct_time): Time stamp of when the optimization
            ended
        message (str): Description of why :func:`.optimize_pulses` completed,
            E.g, "Reached 1000 iterations"

    """

    time_fmt = "%Y-%m-%d %H:%M:%S"
    """Format used in :attr:`start_local_time_str` and
    :attr:`end_local_time_str`
    """

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
        self.message = ''

    def __str__(self):
        return dedent(
            r'''
        Krotov Optimization Result
        --------------------------
        - Started at {start_local_time}
        - Number of objectives: {n_objectives}
        - Number of iterations: {n_iters}
        - Reason for termination: {message}
        - Ended at {end_local_time}
        '''.format(
                start_local_time=self.start_local_time_str,
                n_objectives=len(self.objectives),
                n_iters=len(self.iters) - 1,  # do not count zero iteration
                end_local_time=self.end_local_time_str,
                message=self.message,
            )
        ).strip()

    def __repr__(self):
        return self.__str__()

    @property
    def start_local_time_str(self):
        """The :attr:`start_local_time` attribute formatted as a string"""
        if self.start_local_time is not None:
            return time.strftime(self.time_fmt, self.start_local_time)
        else:
            return 'n/a'

    @property
    def end_local_time_str(self):
        """The :attr:`end_local_time` attribute formatted as a string"""
        if self.end_local_time is not None:
            return time.strftime(self.time_fmt, self.end_local_time)
        else:
            return 'n/a'

    @property
    def optimized_objectives(self):
        """list[Objective]: A copy of the :attr:`objectives` with the
        :attr:`optimized_controls` plugged in"""
        objectives = []
        for (i, obj) in enumerate(self.objectives):
            H = _plug_in_optimized_controls(
                obj.H, self.optimized_controls, self.controls_mapping[i][0]
            )
            c_ops = [
                _plug_in_optimized_controls(
                    c_op,
                    self.optimized_controls,
                    self.controls_mapping[i][j + 1],
                )
                for (j, c_op) in enumerate(obj.c_ops)
            ]
            objectives.append(
                Objective(
                    H=H,
                    initial_state=obj.initial_state,
                    target=obj.target,
                    c_ops=c_ops,
                )
            )
        return objectives

    @classmethod
    def load(cls, filename, objectives=None):
        """Construct :class:`Result` object from a :meth:`dump` file

        Args:
            filename (str): The file from which to load the :class:`Result`.
                Must be in the format created by :meth:`dump`.
            objectives (None or list[Objective]): If given, after loading
                :class:`Result` from the given `filename`, overwrite
                :attr:`objectives` with the given `objectives`. This is
                necessary because :meth:`dump` does not preserve time-dependent
                controls that are Python functions.

        Returns:
            Result: The :class:`Result` instance loaded from `filename`
        """
        with open(filename, 'rb') as dump_fh:
            result = pickle.load(dump_fh)
        if objectives is None:
            for obj in result.objectives:
                if _contains_control_placeholders(obj.H) or any(
                    [_contains_control_placeholders(lst) for lst in obj.c_ops]
                ):
                    logger = logging.getLogger('krotov')
                    logger.warning(
                        "Result.objectives contains control placeholders. "
                        "You should overwrite it."
                    )
                    break
        else:
            result.objectives = objectives
        return result

    def dump(self, filename):
        """Dump the :class:`Result` to a binary :mod:`pickle` file.

        The original :class:`Result` object can be restored from the resulting
        file using :meth:`load`. However, time-dependent control fields that
        are callables/functions will not be preserved, as they are not
        "pickleable".

        Args:
            filename (str): Name of file to which to dump the :class:`Result`.
        """
        with open(filename, 'wb') as dump_fh:
            pickle.dump(self, dump_fh)


def _contains_control_placeholders(lst):
    if isinstance(lst, list):
        return any([_contains_control_placeholders(v) for v in lst])
    else:
        return isinstance(lst, _ControlPlaceholder)


def _plug_in_optimized_controls(H, controls, mapping):
    """Auxilliary routine to :attr:`Result.optimized_objectives`"""
    H = _nested_list_shallow_copy(H)
    for (control, control_mapping) in zip(controls, mapping):
        for i in control_mapping:
            H[i][1] = control
    return H
