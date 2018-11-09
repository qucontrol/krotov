"""Tests for krotov.Result in isolation"""
import numpy as np

import krotov


def test_control_tlist_calculation():
    """Test calculation of control_tlist for non-equidistant time grid"""
    target = krotov.Result([], [], tlist=np.array([0, 1.0, 2.0, 2.2]))
    assert len(target.tlist) == 4
    assert len(target.control_tlist) == len(target.tlist) - 1
    assert target.control_tlist[0] == 0.5
    assert target.control_tlist[1] == 1.5
    assert target.control_tlist[2] == 2.1
