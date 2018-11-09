"""Test conversion of controls between different time grids"""
import numpy as np
from functools import partial

import krotov


def test_conversion_inverse():
    """Test that `controls_onto_interval` and `pulses_onto_tlist` are
    inverses"""
    tlist = np.linspace(0, 10, 20)
    tlist_midpoints = []
    for i in range(len(tlist) - 1):
        tlist_midpoints.append(0.5 * (tlist[i+1] + tlist[i]))
    tlist_midpoints = np.array(tlist_midpoints)

    blackman = partial(krotov.shapes.blackman, t_start=0, t_stop=10)

    pulse_orig = krotov.pulse_conversion.control_onto_interval(
        blackman, tlist, tlist_midpoints)

    control = krotov.pulse_conversion.pulse_onto_tlist(pulse_orig)
    pulse = krotov.pulse_conversion.control_onto_interval(
        control, tlist, tlist_midpoints)

    assert np.max(np.abs(pulse - pulse_orig)) < 1e-14
