"""Routines for structural conversions.

Conversion between between data structures used by QuTiP's
:func:`~qutip.mesolve.mesolve` and data structures used internally in an
optimization with Krotov's method. This includes the
time discretization of control fields, and in particular converting between a
discretization defined on the *points* of the time grid ("controls") and
piecewise-constant "pulses" defined on the *intervals* of the time grid.
"""
import copy
import logging
import warnings

import numpy as np


__all__ = [
    'control_onto_interval',
    'discretize',
    'extract_controls',
    'extract_controls_mapping',
    'plug_in_pulse_values',
    'pulse_onto_tlist',
    'pulse_options_dict_to_list',
]


def _nested_list_shallow_copy(l):
    if isinstance(l, list):
        return [copy.copy(h) if isinstance(h, list) else h for h in l]
    else:
        return l


def _tlist_midpoints(tlist):
    """Calculate array of midpoints in `tlist`."""
    tlist_midpoints = []
    for i in range(len(tlist) - 1):
        tlist_midpoints.append(0.5 * (tlist[i + 1] + tlist[i]))
    return np.array(tlist_midpoints)


def _find_in_list(val, list_to_search):
    """Return index of `val` in `list_to_search`, or -1.

    Works even if `val` is a `numpy.ndarray`. In this case, comparison is by
    object identity.
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


def discretize(control, tlist, args=(None,), kwargs=None, via_midpoints=False):
    """Discretize the given `control` onto the `tlist` time grid.

    If `control` is a callable, return array of values for `control` evaluated
    for all points in `tlist`.  If `control` is already discretized, check that
    the discretization matches `tlist` (by size).

    Args:
        control (callable or numpy.ndarray): control to be discretized. If
            callable, must take time value `t` as its first argument.
        tlist (numpy.ndarray): time grid to discretize one
        args (tuple or list): If `control` is a callable, further positional
            arguments to pass to `control`. The default passes a single value
            None, to match the requirements for a callable control function in
            QuTiP.
        kwargs (None or dict): If `control` is callable, further keyword
            arguments to pass to `control`. If None, no keyword arguments will
            be passed.
        via_midpoints (bool): If True, sample `control` at the midpoints of
            `tlist` (except for the initial and final values which are
            evaluated at ``tlist[0]`` and ``tlist[1]`` to preseve exact
            boundary conditions) and then un-average via
            :func:`pulse_onto_tlist` to fit onto `tlist`. If False, evaluate
            directly on `tlist`.

    Note:
        If ``via_midpoints=True``, the discretized values are generally not
        *exactly* the result of evaluating `control` at the values of `tlist`.
        Instead, the values are adjusted slightly to guarantee numerical
        stability when converting between a sampling on the time grid and a
        sampling on the mid points of the time grid intervals, as required by
        Krotov's method, see :ref:`TimeDiscretization`.

    Returns:
        numpy.ndarray: Discretized array of real `control` values, same length
        as `tlist`

    Raises:
        TypeError: If `control` is not a function that takes two arguments
            (`t`, `args`), or a numpy array
        ValueError: If `control` is numpy array of incorrect size.
    """
    warnings.filterwarnings(action="error", category=np.ComplexWarning)
    # see https://stackoverflow.com/q/54814133/152544
    if callable(control):
        if kwargs is None:
            kwargs = {}
        if via_midpoints:
            tlist_midpoints = (tlist + 0.5 * (tlist[1] - tlist[0]))[:-1]
            tlist_midpoints[0] = tlist[0]
            tlist_midpoints[-1] = tlist[-1]
            pulse_on_midpoints = discretize(
                control,
                tlist_midpoints,
                args=args,
                kwargs=kwargs,
                via_midpoints=False,
            )
            return pulse_onto_tlist(pulse_on_midpoints)
        else:
            # relies on np.ComplexWarning being thrown as an error
            return np.array(
                [float(control(t, *args, **kwargs)) for t in tlist],
                dtype=np.float64,
            )
    elif isinstance(control, (np.ndarray, list)):
        # relies on np.ComplexWarning being thrown as an error
        control = np.array([float(v) for v in control], dtype=np.float64)
        if len(control) != len(tlist):
            raise ValueError(
                "If control is an array, it must of the same length as tlist"
            )
        return control
    else:
        raise TypeError(
            "control must be either a callable func(t, args) or a numpy array"
        )


def extract_controls(objectives):
    """Extract a list of (unique) controls from the `objectives`.

    Controls are unique if they are not the same object, cf.
    `Python's is keyword`_.

    .. _Python's is keyword: https://docs.python.org/3/reference/expressions.html#is

    Args:
        objectives (list): List of :class:`.Objective` instances

    Returns:
        list of controls in `objectives`

    See :func:`extract_controls_mapping` for an example.
    """
    controls = []
    for i_obj, objective in enumerate(objectives):
        for i_ham, ham in enumerate(objective.H):
            if isinstance(ham, list):
                assert len(ham) == 2
                control = ham[1]
                if _find_in_list(control, controls) < 0:
                    controls.append(control)
    return controls


def _control_indices_in_nested_list(nested_list, control):
    """Given a nested list (QuTiP Hamiltonian), find the indices that contain
    `control` and return them as a list"""
    result = []
    for i, item in enumerate(nested_list):
        if isinstance(item, list):
            assert len(item) == 2
            if item[1] is control:
                result.append(i)
    return result


def extract_controls_mapping(objectives, controls):
    """Extract a map of where `controls` are used in `objectives`.

    The result is a nested list where the first index relates to the
    `objectives`, the second index relates to the Hamiltonian (0) or the
    `c_ops` (1...), and the third index relates to the `controls`.

    Example:

        >>> import qutip
        >>> import krotov
        >>> X, Y, Z = qutip.Qobj(), qutip.Qobj(), qutip.Qobj() # dummy Hams
        >>> u1, u2 = np.array([]), np.array([])                # dummy controls
        >>> psi0, psi_tgt = qutip.Qobj(), qutip.Qobj()         # dummy states

        >>> H1 = [X, [Y, u1], [Z, u1]]  # ham for first objective
        >>> H2 = [X, [Y, u2]]           # ham for second objective
        >>> c_ops = [[[X, u1]], [[Y, u2]]]
        >>> objectives = [
        ...     krotov.Objective(
        ...         initial_state=psi0,
        ...         target=psi_tgt,
        ...         H=H1,
        ...         c_ops=c_ops
        ...     ),
        ...     krotov.Objective(
        ...         initial_state=psi0,
        ...         target=psi_tgt,
        ...         H=H2,
        ...         c_ops=c_ops
        ...     )
        ... ]
        >>> controls = extract_controls(objectives)
        >>> assert controls == [u1, u2]

        >>> controls_mapping = extract_controls_mapping(objectives, controls)
        >>> controls_mapping
        [[[[1, 2], []], [[0], []], [[], [0]]], [[[], [1]], [[0], []], [[], [0]]]]

        The structure should be read as follows:

        * For the first objective (0), in the Hamiltonian (0), where is
          the first pulse (0) used? (answer: in ``H1[1]`` and ``H1[2]``)

            >>> controls_mapping[0][0][0]
            [1, 2]

        * For the second objective (1), in the second ``c_ops`` (2), where is
          the second pulse (1) used? (answer: in ``c_ops[1][0]``)

            >>> controls_mapping[1][2][1]
            [0]

        * For the second objective (1), in the Hamiltonian (0), where is the
          first pulse (0) used? (answer: nowhere)

            >>> controls_mapping[1][0][0]
            []
    """
    controls_mapping = []
    for objective in objectives:
        controls_mapping.append([])
        controls_mapping[-1].append(
            [
                _control_indices_in_nested_list(objective.H, control)
                for control in controls
            ]
        )
        for c_op in objective.c_ops:
            controls_mapping[-1].append(
                [
                    _control_indices_in_nested_list(c_op, control)
                    for control in controls
                ]
            )
    return controls_mapping


def pulse_options_dict_to_list(pulse_options, controls):
    """Convert `pulse_options` into a list.

    Given a dict `pulse_options` that contains an options-dict
    for every control in `controls` (cf. :func:`.optimize_pulses`), return a
    list of the options-dicts in the same order as `controls`.

    Raises:
        ValueError: if `pulse_options` to not contain all of the `controls`
    """
    logger = logging.getLogger('krotov')
    if len(pulse_options) > len(controls):
        logger.warning(
            "pulse_options contains extra elements that are not in `controls`"
        )
    pulse_options_list = []
    for control in controls:
        try:
            try:
                opts = pulse_options[control]
            except TypeError:  # control is numpy array
                opts = pulse_options[id(control)]
            pulse_options_list.append(opts)
        except KeyError:
            raise ValueError(
                "The control %s does not have any associated pulse options"
                % str(control)
            )
    return pulse_options_list


def plug_in_pulse_values(H, pulses, mapping, time_index, conjugate=False):
    """Plug pulse values into H.

    Args:
        H (list): nested list for a QuTiP-time-dependent operator
        pulses (list): list of pulses in array format
        mapping (list): nested list: for each pulse, a list of indices in `H`
            where pulse value should be inserted
        time_index (int): Index of the value of each pulse that should be
            plugged in
        conjugate (bool): If True, use conjugate complex pulse values

    Returns:
        list: a list with the same structure as `H` that contains the same
        :class:`~qutip.Qobj` operators as `H`, but where every time dependency
        is replaced by the value of the appropriate pulse at `time_index`.

    Example:

        >>> X, Y, Z = 'X', 'Y', 'Z' # dummy Hams, these would normally be Qobjs
        >>> u1, u2 = np.array([0, 10, 0]), np.array([0, 20, 0])
        >>> H = [X, [X, u1], [Y, u1], [Z, u2]]
        >>> pulses = [u1, u2]
        >>> mapping = [[1, 2], [3]]  # u1 is in H[1] and H[2], u2 is in H[3]
        >>> plug_in_pulse_values(H, pulses, mapping, time_index=1)
        ['X', ['X', 10], ['Y', 10], ['Z', 20]]

    Note:
        It is of no consequence whether `H` contains the `pulses`, as long as
        it has the right structure::

            >>> H = [X, [X, None], [Y, None], [Z, None]]
            >>> plug_in_pulse_values(H, pulses, mapping, time_index=1)
            ['X', ['X', 10], ['Y', 10], ['Z', 20]]
    """
    H = _nested_list_shallow_copy(H)
    for (pulse, pulse_mapping) in zip(pulses, mapping):
        for i in pulse_mapping:
            if conjugate:
                H[i][1] = np.conjugate(pulse[time_index])
            else:
                H[i][1] = pulse[time_index]
    return H


def control_onto_interval(control):
    """Convert control on time grid to control on time grid intervals.

    Args:
        control (numpy.ndarray): values of controls on time grid

    Returns:
        numpy.ndarray: pulse defined on the intervals of the time grid

    The value for the first and last interval will be identical to the values
    at ``control[0]`` and ``control[-1]`` to ensure proper boundary conditions.
    All other intervals are calculated such that the original values in
    `control` are the average of the interval-values before and after that
    point in time.

    The :func:`pulse_onto_tlist` function calculates the inverse to this
    transformation.

    Note:
        For a callable `control`, call :func:`discretize` first.
    """
    if isinstance(control, np.ndarray):
        assert len(control.shape) == 1  # must be 1D array
        pulse = np.zeros(len(control) - 1, dtype=control.dtype.type)
        pulse[0] = control[0]
        for i in range(1, len(control) - 1):
            pulse[i] = 2.0 * control[i] - pulse[i - 1]
        pulse[-1] = control[-1]
        return pulse
    else:
        raise ValueError(
            "Not implemented: control type %s" % control.__class__.__name__
        )


def pulse_onto_tlist(pulse):
    """Convert `pulse` from time-grid intervals to time-grid points.

    Args:
        pulse (numpy.ndarray): values defined on the interval of a time grid

    Returns:
        numpy.ndarray: values of the control defined directly on the time grid
        points. The size of the returned array is one greater than the size
        of `pulse`.

    Inverse of :func:`control_onto_interval`.

    The first and last value are also the first and last value of the returned
    control field. For all other points, the value is the average of the value
    of the input values before and after the point.
    """
    control = np.zeros(len(pulse) + 1, dtype=pulse.dtype.type)
    control[0] = pulse[0]
    for i in range(1, len(control) - 1):
        control[i] = 0.5 * (pulse[i - 1] + pulse[i])
    control[-1] = pulse[-1]
    return control
