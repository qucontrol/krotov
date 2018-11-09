"""Tests for krotov.Result in isolation"""
import numpy as np

import krotov


def test_control_tlist_calculation():
    """Test calculation of control_tlist for non-equidistant time grid"""
    target = krotov.Result([], [], tlist=np.array([0, 1.0, 2.0, 2.2]))
    assert len(target.tlist) == 4
    assert len(target.tlist_midpoints) == len(target.tlist) - 1
    assert target.tlist_midpoints[0] == 0.5
    assert target.tlist_midpoints[1] == 1.5
    assert target.tlist_midpoints[2] == 2.1
