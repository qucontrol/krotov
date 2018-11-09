"""Functions that may be used for the `shape` attribute of
:class:`.PulseOptions`, or for generating guess pulses"""

import numpy as np

__all__ = ['flattop', 'box', 'blackman']


def flattop(t, t_start, t_stop, t_rise, t_fall=None, func='blackman'):
    """Flat shape (one) with a switch-on/switch-off from zero

    The flattop function starts at 0, and ramps to to 1 during the `t_rise`
    interval. For ``func='blackman'``, the switch-on shape is half of a
    Blackman window (see :func:`blackman`). For ``func='sinsq``, it is a
    sine-squared curve. The function then remains at value 1, before ramping
    down to 0 again during `t_fall`.

    Args:
        t (float): Time  point or time grid
        t_start (float): Start of flattop window
        t_stop (float): Stop of flattop window
        t_rise (float): Duration of ramp-up, starting at `t_start`
        t_fall (float): Duration of ramp-down, ending at `t_stop`.
            If not given, ``t_fall=t_rise``.
        func (str): One of 'blackman', 'sinsq'

    Note:
        You may use :class:`numpy.vectorize` to transform this into a shape
        function for arrays, and :func:`functools.partial` to fix the function
        arguments other than `t`, creating a function suitable for the `shape`
        attribute of :class:`.PulseOptions`.
    """
    if t_fall is None:
        t_fall = t_rise
    if func == 'blackman':
        return _flattop_blackman(t, t_start, t_stop, t_rise, t_fall)
    elif func == 'sinsq':
        return _flattop_sinsq(t, t_start, t_stop, t_rise, t_fall)
    else:
        raise ValueError("Invalid func: %s" % func)


def _flattop_sinsq(t, t_start, t_stop, t_rise, t_fall):
    if (t >= t_start) and (t <= t_stop):
        f = 1.0
        if t <= t_start + t_rise:
            f = np.sin(np.pi * (t-t_start) / (2.0*t_rise))**2
        elif t >= t_stop - t_fall:
            f = np.sin(np.pi * (t-t_stop) / (2.0*t_fall))**2
        return f
    else:
        return 0.0


def _flattop_blackman(t, t_start, t_stop, t_rise, t_fall):
    if (t >= t_start) and (t <= t_stop):
        f = 1.0
        if t <= t_start + t_rise:
            f = blackman(t, t_start, t_start + 2*t_rise)
        elif t >= t_stop - t_fall:
            f = blackman(t, t_stop - 2*t_fall, t_stop)
        return f
    else:
        return 0.0


def box(t, t_start, t_stop):
    """Box-shape (Theta-function)

    The shape is 0 before `t_start` and after `t_stop` and 1 elsewhere.

    Args:
        t (float): Time point or time grid
        t_start (float): First value of `t` for which the box has value 1
        t_stop (float): Last value of `t` for which the box has value 1

    Note:
        You may use :class:`numpy.vectorize` and :func:`functools.partial`, cf.
        :func:`flattop`.
    """
    if t < t_start:
        return 0.0
    if t > t_stop:
        return 0.0
    return 1.0


def blackman(t, t_start, t_stop, a=0.16):
    """Blackman window shape

    See http://en.wikipedia.org/wiki/Window_function#Blackman_windows

    A Blackman shape looks nearly identical to a Gaussian with a 6-sigma
    interval between start and stop  Unlike the Gaussian,
    however, it will go exactly to zero at the edges. Thus, Blackman pulses
    are often preferable to Gaussians.

    Args:
        t (float or numpy array): Time point or time grid
        t_start (float): Starting point of Blackman shape
        t_stop (float): End point of Blackman shape

    Returns:
        float or numpy array: If `t` is a float, return the value of the
        Blackman shape at `t`.  If `t` is an array, return an array of same
        size as `t`, containing the values for the Blackman shape (zero before
        `t_start` and after `t_stop`)
    """
    T = t_stop - t_start
    box_vec = np.vectorize(box)
    return (
        0.5 * box_vec(t, t_start, t_stop) * (
            1.0 - a - np.cos(2.0*np.pi * (t-t_start)/T) +
            a * np.cos(4.0 * np.pi * (t-t_start)/T)))
