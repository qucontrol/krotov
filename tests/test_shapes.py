"""Test of krotov.shape functions"""
from functools import partial

import pytest

from krotov.shapes import flattop


def test_flattop_blackman():
    """Check basic properties of a Blackman flattop shape"""
    shape = partial(flattop, t_start=10, t_stop=20, t_rise=2, func='blackman')
    assert shape(9.9) == 0
    assert shape(10) < 1e-14
    assert shape(20) < 1e-14
    assert shape(20.1) == 0
    assert shape(15) == 1


def test_flattop_sinsq():
    """Check basic properties of a sinsq flattop shape"""
    shape = partial(flattop, t_start=10, t_stop=20, t_rise=2, func='sinsq')
    assert shape(9.9) == 0
    assert shape(10) < 1e-14
    assert shape(20) < 1e-14
    assert shape(20.1) == 0
    assert shape(15) == 1


def test_invalid_flattop():
    """Check that a flattop with an invalid 'func' raises an Exception"""
    with pytest.raises(ValueError):
        flattop(0, t_start=10, t_stop=20, t_rise=2, func='xxx')
