"""Test of structural conversions"""
import logging
from functools import partial

import numpy as np
import pytest
import qutip

import krotov
from krotov.conversions import (
    discretize,
    extract_controls,
    extract_controls_mapping,
    pulse_options_dict_to_list,
)
from krotov.shapes import qutip_callback


def test_conversion_control_pulse_inverse():
    """Test that `controls_onto_interval` and `pulses_onto_tlist` are
    inverses"""
    tlist = np.linspace(0, 10, 20)

    blackman = qutip_callback(krotov.shapes.blackman, t_start=0, t_stop=10)

    pulse_orig = krotov.conversions.control_onto_interval(
        discretize(blackman, tlist)
    )

    control = krotov.conversions.pulse_onto_tlist(pulse_orig)
    pulse = krotov.conversions.control_onto_interval(control)

    assert np.max(np.abs(pulse - pulse_orig)) < 1e-14


def test_discretize():
    """Test the `discretize` routine"""
    tlist = np.linspace(0, 10, 20)

    control = partial(krotov.shapes.blackman, t_start=0, t_stop=10)
    with pytest.raises(TypeError):
        discretize(control, tlist)

    with pytest.raises(TypeError):
        discretize('sin(t)', tlist)

    control = qutip_callback(krotov.shapes.blackman, t_start=0, t_stop=10)

    with pytest.raises(ValueError):
        discretize(np.array([control(t, None) for t in tlist[:-1]]), tlist)

    control_array = discretize(control, tlist)
    assert len(control_array) == len(tlist)
    assert abs(control_array[0]) < 1e-15
    assert abs(control_array[-1]) < 1e-15

    control = np.array([control(t, None) for t in tlist])
    assert np.max(np.abs((control - control_array))) < 1e-15

    control_array2 = discretize(control, tlist)
    assert np.max(np.abs((control_array2 - control_array))) < 1e-15


def test_discretization_as_float():
    """Test that a control function that returns an int still gets discretized
    as a float (#41)"""
    H = [qutip.Qobj(), [qutip.Qobj(), lambda t, args: 0]]
    u = H[1][1]
    objectives = [
        krotov.Objective(initial_state=qutip.Qobj(), target=None, H=H)
    ]
    tlist = np.linspace(0, 10, 100)

    res = krotov.optimize._initialize_krotov_controls(
        objectives, {u: dict(lambda_a=1, update_shape=lambda t: 0)}, tlist
    )

    guess_controls = res[0]
    assert guess_controls[0].dtype == np.float64
    guess_pulses = res[1]
    assert guess_pulses[0].dtype == np.float64
    shape_arrays = res[4]
    assert shape_arrays[0].dtype == np.float64


def test_initialize_krotov_controls():
    """Check that pulses and controls are initialized while preserving the
    correct boundary conditions.

    This is the point that the section "Time Discretization Schemes" in the
    documentation is making.

    Tests the resolution of #20.
    """

    T = 10
    blackman = qutip_callback(krotov.shapes.blackman, t_start=0, t_stop=T)
    H = ['H0', ['H1', blackman]]
    tlist = np.linspace(0, T, 10)
    pulse_options = {blackman: dict(lambda_a=1.0, update_shape=1)}

    objectives = [
        krotov.Objective(initial_state=qutip.Qobj(), target=None, H=H)
    ]

    assert abs(blackman(0, None)) < 1e-15
    assert abs(blackman(T, None)) < 1e-15

    (
        guess_controls,
        guess_pulses,
        pulses_mapping,
        lambda_vals,
        shape_arrays,
    ) = krotov.optimize._initialize_krotov_controls(
        objectives, pulse_options, tlist
    )

    assert isinstance(guess_controls[0], np.ndarray)
    assert len(guess_controls[0]) == len(tlist)
    assert abs(guess_controls[0][0]) < 1e-15
    assert abs(guess_controls[0][-1]) < 1e-15

    assert isinstance(guess_pulses[0], np.ndarray)
    assert len(guess_pulses[0]) == len(tlist) - 1
    assert abs(guess_pulses[0][0]) < 1e-15
    assert abs(guess_pulses[0][-1]) < 1e-15

    assert len(pulse_options) == 1

    assert len(pulses_mapping) == 1
    assert len(pulses_mapping[0]) == 1
    assert len(pulses_mapping[0][0]) == 1
    assert len(pulses_mapping[0][0][0]) == 1
    assert pulses_mapping[0][0][0][0] == 1

    assert len(lambda_vals) == 1
    assert lambda_vals[0] == 1.0

    assert len(shape_arrays) == 1
    assert isinstance(shape_arrays[0], np.ndarray)
    assert len(shape_arrays[0]) == len(tlist) - 1


def test_extract_controls_with_arrays():
    """Test extract_controls for controls that are numpy arrays"""
    X, Y, Z = qutip.Qobj(), qutip.Qobj(), qutip.Qobj()  # dummy Hamiltonians
    u1, u2 = np.array([]), np.array([])  # dummy controls
    psi0, psi_tgt = qutip.Qobj(), qutip.Qobj()  # dummy states

    assert X is not Y
    assert Y is not Z
    assert u1 is not u2
    assert psi0 is not psi_tgt

    H1 = [X, [Y, u1], [Z, u2]]  # ham for first objective
    H2 = [X, [Y, u2]]  # ham for second objective
    objectives = [
        krotov.Objective(initial_state=psi0, target=psi_tgt, H=H1),
        krotov.Objective(initial_state=psi0, target=psi_tgt, H=H2),
    ]

    controls = extract_controls(objectives)
    control_map = extract_controls_mapping(objectives, controls)

    assert controls == [u1, u2]
    assert control_map[0] == [[[1], [2]]]
    assert control_map[1] == [[[], [1]]]


def test_extract_controls():
    """Check that we can extract a list of controls from a list of
    objectives"""

    # dummy objects
    X = qutip.Qobj()
    Y = qutip.Qobj()
    f = lambda t: 0
    g = lambda t: 0
    h = lambda t: 0
    d = lambda t: 0

    # all of these dummy objects should be distinguishable
    assert f is not g
    assert g is not h
    assert h is not d
    assert X is not Y

    H1 = [X, [X, f], [X, g]]
    H2 = [X, [X, f], [X, h]]
    H3 = [X, [X, d], X]

    # check same Hamiltonian occuring in multiple objectives
    objectives = [
        krotov.Objective(initial_state=X, target=Y, H=H1),
        krotov.Objective(initial_state=Y, target=X, H=H1),
    ]
    controls = extract_controls(objectives)
    maps = extract_controls_mapping(objectives, controls)
    assert len(controls) == 2
    assert f in controls
    assert g in controls
    assert len(maps) == len(objectives)
    assert maps == [[[[1], [2]]], [[[1], [2]]]]

    # check same control occuring in multiple Hamiltonians
    objectives = [
        krotov.Objective(initial_state=X, target=Y, H=H1),
        krotov.Objective(initial_state=Y, target=X, H=H2),
        krotov.Objective(initial_state=Y, target=X, H=H3),
    ]
    controls = extract_controls(objectives)
    assert len(controls) == 4
    for c in (f, g, h, d):
        assert c in controls
    maps = extract_controls_mapping(objectives, controls)
    assert maps[0] == [[[1], [2], [], []]]
    assert maps[1] == [[[1], [], [2], []]]
    assert maps[2] == [[[], [], [], [1]]]


def test_pulse_options_dict_to_list(caplog):
    """Test conversion of pulse_options"""

    u1, u2, u3 = np.array([]), np.array([]), np.array([])  # dummy controls
    controls = [u1, u2]

    assert u1 is not u2
    assert u2 is not u3

    pulse_options = {
        id(u1): dict(lambda_a=1.0, update_shape=1),
        id(u2): dict(lambda_a=2.0, update_shape=1),
    }

    pulse_options_list = pulse_options_dict_to_list(pulse_options, controls)
    assert len(pulse_options_list) == 2
    assert pulse_options_list[0] == pulse_options[id(u1)]
    assert pulse_options_list[1] == pulse_options[id(u2)]

    # check error for missing pulse options
    pulse_options = {id(u1): dict(lambda_a=1.0, update_shape=1)}
    with pytest.raises(ValueError) as exc_info:
        pulse_options_dict_to_list(pulse_options, controls)
    assert 'does not have any associated pulse options' in str(exc_info.value)

    # check warning message for extra pulse options
    pulse_options = {
        id(u1): dict(lambda_a=1.0, update_shape=1),
        id(u2): dict(lambda_a=1.0, update_shape=1),
        id(u3): dict(lambda_a=1.0, update_shape=1),
    }
    with caplog.at_level(logging.WARNING):
        pulse_options_dict_to_list(pulse_options, controls)
    assert 'extra elements' in caplog.text


def test_control_tlist_calculation():
    """Test calculation of tlist_midpoints for non-equidistant time grid"""
    tlist = np.array([0, 1.0, 2.0, 2.2])
    midpoints = krotov.conversions._tlist_midpoints(tlist)
    assert len(midpoints) == len(tlist) - 1
    assert midpoints[0] == 0.5
    assert midpoints[1] == 1.5
    assert midpoints[2] == 2.1
