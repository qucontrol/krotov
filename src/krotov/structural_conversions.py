"""Routines for converting converting between structures good for QuTiP's
mesolve and Krotov"""
import numpy as np
import logging

__all__ = [
    'control_onto_interval', 'pulse_onto_tlist', 'extract_controls',
    'extract_controls_mapping', 'pulse_options_dict_to_list']


def _tlist_midpoints(tlist):
    """Calculate array of midpoints in `tlist`"""
    tlist_midpoints = []
    for i in range(len(tlist) - 1):
        tlist_midpoints.append(0.5 * (tlist[i+1] + tlist[i]))
    return np.array(tlist_midpoints)


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


def extract_controls(objectives):
    """Extract a list of (unique) controls from the `objectives`

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
    """Extract a map of where `controls` are used in `objectives`

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
        >>> c_ops = ([[X, u1]], [[Y, u2]])
        >>> objectives = [
        ...     krotov.Objective(H1, psi0, psi_tgt, c_ops=c_ops),
        ...     krotov.Objective(H2, psi0, psi_tgt, c_ops=c_ops)]
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
        controls_mapping[-1].append([
            _control_indices_in_nested_list(objective.H, control)
            for control in controls])
        for c_op in objective.c_ops:
            controls_mapping[-1].append([
                _control_indices_in_nested_list(c_op, control)
                for control in controls])
    return controls_mapping


def pulse_options_dict_to_list(pulse_options, controls):
    """Convert `pulse_options` into a list

    Given a dict `pulse_options` that contains a :class:`.PulseOptions`
    instance for every control in `controls`, return a list of the
    :class:`.PulseOptions` in the same order as `controls`.

    Raises:
        ValueError: if `pulse_options` to not contain all of the `controls`
    """
    logger = logging.getLogger('krotov')
    if len(pulse_options) > len(controls):
        logger.warning(
            "pulse_options contains extra elements that are not in `controls`")
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
                % str(control))
    return pulse_options_list


def control_onto_interval(control, tlist, tlist_midpoints):
    """Convert control on `tlist` to `tlist` intervals (`tlist_midpoints`)

    Args:
        control (callable or numpy array): values at `tlist`, either as a
            function ``control(t, args)`` or an array of the same length as
            `tlist`
        tlist (numpy array): time grid point values
        tlist_midpoints (numpy array): midpoint values in `tlist_midpoints`.

    Returns:
        numpy array: pulse defined on the intervals to `tlist`, that is the
        `tlist_midpoints`.

    The value for the first and last interval will be identical to the values
    at ``tlist[0]`` and ``tlist[-1]`` to ensure proper boundary conditions. All
    other intervals are calculated such that the original values on `tlist` are
    the average of the interval-values before and after that point in time.

    The :func:`pulse_onto_tlist` function calculates the inverse to this
    transformation.
    """
    assert len(tlist_midpoints) == len(tlist) - 1
    if callable(control):
        return np.array([control(t, None) for t in tlist_midpoints])
    elif isinstance(control, np.ndarray):
        assert len(control) == len(tlist)
        assert len(control.shape) == 1  # must be 1D array
        pulse = np.zeros(len(control)-1, dtype=control.dtype.type)
        pulse[0] = control[0]
        for i in range(1, len(control)-1):
            pulse[i] = 2.0 * control[i] - pulse[i-1]
        pulse[-1] = control[-1]
        return pulse
    else:
        raise ValueError(
            "Not implemented: control type %s"
            % control.__class__.__name__)


def pulse_onto_tlist(pulse):
    """Convert `pulse` from time-grid intervals to time-grid points

    Args:
        pulse (numpy array): values defined on the interval of a time grid

    Returns:
        numpy array: values of the control defined directly on the time grid
        points. The size of the returned array is one greater than the size
        of `pulse`.

    Inverse of :func:`control_onto_interval`.

    The first and last value are also the first and last value of the returned
    control field. For all other points, the value is the average of the value
    of the input values before and after the point.
    """
    control = np.zeros(len(pulse)+1, dtype=pulse.dtype.type)
    control[0] = pulse[0]
    for i in range(1, len(control)-1):
        control[i] = 0.5 * (pulse[i-1] + pulse[i])
    control[-1] = pulse[-1]
    return control
