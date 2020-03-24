"""Module defining the :class:`Result` object that is returned by
:func:`.optimize_pulses`.
"""
import copyreg
import datetime
import logging
import pickle
import time
from textwrap import dedent

from .conversions import _nested_list_shallow_copy, pulse_onto_tlist
from .objectives import Objective, _ControlPlaceholder, _Objective_reduce


__all__ = ['Result']


class Result:
    """Result of a Krotov optimization with :func:`.optimize_pulses`.

    Attributes:
        objectives (list[Objective]): The control objectives
        tlist (numpy.ndarray): The time grid values
        iters (list[int]): Iteration numbers, starting at 0.
        iter_seconds (list[int]): for each iteration number, the number of
            seconds that were spent in the optimization
        info_vals (list): For each iteration, the return value of `info_hook`,
            or None
        tau_vals (list[list[complex]): for each iteration, a list of complex
            overlaps between the target state and the forward-propagated state
            for each objective, assuming :attr:`.Objective.target` contains the
            target state. If there is no target state, an empty list.
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
        - Ended at {end_local_time} ({time_delta})
        '''.format(
                start_local_time=self.start_local_time_str,
                n_objectives=len(self.objectives),
                n_iters=len(self.iters) - 1,  # do not count zero iteration
                end_local_time=self.end_local_time_str,
                time_delta=str(
                    datetime.timedelta(
                        seconds=time.mktime(self.end_local_time)
                        - time.mktime(self.start_local_time)
                    )
                ),
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
        :attr:`optimized_controls` plugged in."""
        return self.objectives_with_controls(self.optimized_controls)

    def objectives_with_controls(self, controls):
        """List of objectives with the given `controls` plugged in.

        Args:
            controls (list[numpy.ndarray]): A list of control fields, defined
                on the points of :attr:`tlist`. Must be of the same length as
                :attr:`guess_controls` and :attr:`optimized_controls`.

        Returns:
            list[Objective]: A copy of :attr:`objectives`, where all control
            fields are replaced by the elements of the `controls`.

        Raises:
            ValueError: If `controls` does not have the same number controls as
                :attr:`guess_controls` and :attr:`optimized_controls`, or if
                any `controls` are not defined on the points of the time grid.

        See also:
            For plugging in the optimized controls, the
            :attr:`optimized_objectives` attribute is equivalent to
            ``result.objectives_with_controls(result.optimized_controls)``.
        """
        n = len(self.guess_controls)
        m = len(controls)
        if n != m:
            raise ValueError("Expected %d controls, %d given" % (n, m))
        for control in controls:
            try:
                if len(control) != len(self.tlist):
                    raise ValueError(
                        "controls are not defined on the points of the time "
                        "grid: control has %d values for %d time grid points"
                        % (len(control), len(self.tlist))
                    )
            except TypeError:
                pass  # control is not a numpy array (but maybe a callable)
        objectives = []
        for (i, obj) in enumerate(self.objectives):
            H = _plug_in_optimized_controls(
                obj.H, controls, self.controls_mapping[i][0]
            )
            c_ops = [
                _plug_in_optimized_controls(
                    c_op, controls, self.controls_mapping[i][j + 1]
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
    def load(cls, filename, objectives=None, finalize=False):
        """Construct :class:`Result` object from a :meth:`dump` file

        Args:
            filename (str): The file from which to load the :class:`Result`.
                Must be in the format created by :meth:`dump`.
            objectives (None or list[Objective]): If given, after loading
                :class:`Result` from the given `filename`, overwrite
                :attr:`objectives` with the given `objectives`. This is
                necessary because :meth:`dump` does not preserve time-dependent
                controls that are Python functions.
            finalize (bool): If given as True, make sure that the
                :attr:`optimized_controls` are properly finalized. This allows
                to load a :class:`Result` that was dumped before
                :func:`.optimize_pulses` finished, e.g. by
                :func:`.dump_result`.

        Returns:
            Result: The :class:`Result` instance loaded from `filename`
        """
        logger = logging.getLogger('krotov')
        with open(filename, 'rb') as dump_fh:
            result = pickle.load(dump_fh)
        if objectives is None:
            for obj in result.objectives:
                if _contains_control_placeholders(obj.H) or any(
                    [_contains_control_placeholders(lst) for lst in obj.c_ops]
                ):
                    logger.warning(
                        "Result.objectives contains control placeholders. "
                        "You should overwrite it by passing `objectives`."
                    )
                    break  # only warn once
        else:
            result.objectives = objectives
        nt = len(result.tlist)
        for i, control in enumerate(result.optimized_controls):
            if len(control) == nt:
                pass  # all ok
            elif len(control) == nt - 1:
                if finalize:
                    result.optimized_controls[i] = pulse_onto_tlist(control)
                else:
                    logger.warning(
                        "Result.optimized_controls are not finalized. "
                        "Consider loading with `finalize=True`."
                    )
                    break  # only warn once
            else:
                logger.error(
                    "Result.optimized_controls are incongruent with "
                    "Result.tlist"
                )
                break  # only warn once
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
            pickler = pickle.Pickler(dump_fh)
            pickler.dispatch_table = copyreg.dispatch_table.copy()
            pickler.dispatch_table[Objective] = _Objective_reduce
            pickler.dump(self)


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
