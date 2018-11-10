"""Routines for converting converting between structures good for QuTiP's
mesolve and Krotov"""
import numpy as np
import logging
import copy

__all__ = ['control_onto_interval', 'pulse_onto_tlist']


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


def extract_controls(objectives, pulse_options):
    """Extract a list of controls from the `objectives`, modify the
    Hamiltonians into an "extracted" form, and check that all controls have
    proper options

    Args:
        objectives (list): List of :class:`Objective` instances
        pulse_options (dict): Mapping of time-dependent controls found in the
            Hamiltonians of the objectives to :class:`PulseOptions` instances.
            There must be a mapping for each control.

    Returns:
        tuple: A tuple of the following:

        - list of `controls` extracted from `objective`
        - list of "pointers" for where each control occurs in the `objectives`.
          Each element is a list of tuples (objective-index,
          ham-component-index), where the objective-index gives the index of an
          objective that contains the control, and ham-component-index gives
          the index of a component of the Hamiltonian that is linear in the
          control
        - list of values from `pulse_options`, in the same order as `controls`
        - modified copy of `objectives`, where for each Hamiltonian component
          that contains a control, that control has been replaced with an
          integer that specifies in index of the corresponding control in
          `controls`.

    Example:

        >>> import qutip
        >>> import krotov
        >>> X, Y, Z = qutip.Qobj(), qutip.Qobj(), qutip.Qobj() # dummy Hams
        >>> u1, u2 = np.array([]), np.array([])                # dummy controls
        >>> psi0, psi_tgt = qutip.Qobj(), qutip.Qobj()         # dummy states

        >>> H1 = [X, [Y, u1], [Z, u2]]  # ham for first objective
        >>> H2 = [X, [Y, u2]]           # ham for second objective
        >>> objectives = [
        ...     krotov.Objective(H1, psi0, psi_tgt),
        ...     krotov.Objective(H2, psi0, psi_tgt)]
        >>> pulse_options = {
        ...     id(u1): krotov.PulseOptions(lambda_a=1.0),
        ...     id(u2): krotov.PulseOptions(lambda_a=1.0)}

        >>> controls, control_map, options, objectives = extract_controls(
        ...     objectives, pulse_options)
        >>> assert controls == [u1, u2]
        >>> assert objectives[0].H == [X, [Y, 0], [Z, 1]]
        >>> assert objectives[1].H == [X, [Y, 1]]
        >>> assert control_map[0] == [(0, 1)]          # where does u1 occur?
        >>> assert control_map[1] == [(0, 2), (1, 1)]  # where does u2 occur?
        >>> assert options[0] == pulse_options[id(u1)]
        >>> assert options[1] == pulse_options[id(u2)]
    """
    logger = logging.getLogger('krotov')
    controls = []
    controls_map = []
    options_list = []
    objectives = [copy.copy(o) for o in objectives]
    for i_obj, objective in enumerate(objectives):
        for i_ham, ham in enumerate(objective.H):
            if isinstance(ham, list):
                assert len(ham) == 2
                control = ham[1]
                i_control = _find_in_list(control, controls)
                if i_control >= 0:
                    controls_map[i_control].append((i_obj, i_ham))
                else:  # this is a control we haven't seen before
                    try:
                        try:
                            options_list.append(pulse_options[control])
                        except TypeError:  # control is numpy array
                            options_list.append(pulse_options[id(control)])
                    except KeyError:
                        raise ValueError(
                            "The control %s in the component %d of the "
                            "Hamiltonian of the objective %d does not have "
                            "any associated pulse options"
                            % (control, i_ham, i_obj))
                    controls.append(control)
                    controls_map.append([(i_obj, i_ham)])
                    i_control = len(controls) - 1
                ham[1] = i_control
    if len(controls) != len(pulse_options):
        logger.warning(
            "pulse_options contains options for controls that are not in the "
            "objectives")
    return controls, controls_map, options_list, objectives


def control_onto_interval(control, tlist, tlist_midpoints):
    """Convert control on `tlist` to `tlist` intervals (`tlist_midpoints`)

    Args:
        control (callable or numpy array): values at `tlist`, either as a
            function ``control(t, args)`` or an an array of the same length as
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
