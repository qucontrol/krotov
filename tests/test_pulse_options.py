"""Tests of PulseOptions"""
import numpy as np
import pytest
import qutip

import krotov
from krotov.optimize import _initialize_krotov_controls


def test_shape_validation():
    """Test that OCT pulse shapes are converted and verified correctly"""

    H = [qutip.Qobj(), [qutip.Qobj(), lambda t, args: 0]]
    u = H[1][1]
    objectives = [
        krotov.Objective(initial_state=qutip.Qobj(), target=None, H=H)
    ]
    tlist = np.linspace(0, 10, 100)

    res = _initialize_krotov_controls(
        objectives, {u: dict(lambda_a=1, update_shape=1)}, tlist
    )
    # res consists of:
    # guess_controls, guess_pulses, pulses_mapping, lambda_vals, shape_arrays
    shape_arrays = res[4]
    shape_array = shape_arrays[0]
    assert len(shape_arrays) == 1
    assert len(shape_array) == len(tlist) - 1
    assert np.all(shape_array == 1)
    lambda_vals = res[3]
    assert len(lambda_vals) == 1
    assert lambda_vals[0] == 1
    assert isinstance(lambda_vals[0], float)

    res = _initialize_krotov_controls(
        objectives, {u: dict(lambda_a=1, update_shape=0)}, tlist
    )
    shape_array = res[4][0]
    assert np.all(shape_array == 0)

    with pytest.raises(ValueError) as exc_info:
        _initialize_krotov_controls(objectives, {u: dict(lambda_a=1)}, tlist)
    assert "key 'update_shape'" in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        _initialize_krotov_controls(
            objectives, {u: {'update_shape': 1}}, tlist
        )
    assert "key 'lambda_a'" in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        _initialize_krotov_controls(
            objectives, {u: dict(lambda_a=1, update_shape=2)}, tlist
        )
    assert 'update_shape must be a callable' in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        _initialize_krotov_controls(
            objectives,
            {u: dict(lambda_a=1, update_shape=lambda t: 2.0)},
            tlist,
        )
    assert 'in the range [0, 1]' in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        _initialize_krotov_controls(
            objectives,
            {u: dict(lambda_a=1, update_shape=lambda t: 0.5j)},
            tlist,
        )
    assert 'real-valued' in str(exc_info.value)
