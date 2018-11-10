"""Tests of PulseOptions"""
import numpy as np
import pytest

import krotov


def test_shape_validation():
    """Test that OCT pulse shapes are converted and verified correctly"""
    opt = krotov.PulseOptions(lambda_a=1)
    assert callable(opt.shape)
    assert opt.shape(0) == 1

    opt = krotov.PulseOptions(lambda_a=1, shape=0)
    assert callable(opt.shape)
    assert opt.shape(0) == 0

    opt = krotov.PulseOptions(lambda_a=1, shape=1)
    assert callable(opt.shape)
    assert opt.shape(0) == 1

    with pytest.raises(ValueError):
        krotov.PulseOptions(lambda_a=1, shape=2)

    with pytest.raises(ValueError):
        krotov.PulseOptions(lambda_a=1, shape=np.array([0, 0.5, 1, 0.5, 0]))
