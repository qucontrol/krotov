"""Routines for converting pulses between time-grid intervals and time-grid
points"""
import numpy as np


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
