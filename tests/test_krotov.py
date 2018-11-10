"""High-level tests for `krotov` package."""
import qutip
import numpy as np
import logging

import pytest
from pkg_resources import parse_version

import krotov


def test_valid_version():
    """Check that the package defines a valid __version__"""
    assert parse_version(krotov.__version__) >= parse_version("0.0.1")


def test_extract_controls_with_arrays():
    """Test _extract_controls for controls that are numpy arrays"""
    X, Y, Z = qutip.Qobj(), qutip.Qobj(), qutip.Qobj()  # dummy Hamiltonians
    u1, u2 = np.array([]), np.array([])                 # dummy controls
    psi0, psi_tgt = qutip.Qobj(), qutip.Qobj()          # dummy states

    assert X is not Y
    assert Y is not Z
    assert u1 is not u2
    assert psi0 is not psi_tgt

    H1 = [X, [Y, u1], [Z, u2]]  # ham for first objective
    H2 = [X, [Y, u2]]           # ham for second objective
    objectives = [
        krotov.Objective(H1, psi0, psi_tgt),
        krotov.Objective(H2, psi0, psi_tgt)]
    pulse_options = {
        id(u1): krotov.PulseOptions(lambda_a=1.0),
        id(u2): krotov.PulseOptions(lambda_a=1.0)}

    controls, control_map, opts, objectives = krotov._extract_controls(
        objectives, pulse_options)
    assert controls == [u1, u2]
    assert objectives[0].H == [X, [Y, 0], [Z, 1]]
    assert objectives[1].H == [X, [Y, 1]]
    assert control_map[0] == [(0, 1)]          # where does u1 occur?
    assert control_map[1] == [(0, 2), (1, 1)]  # where does u2 occur?
    assert len(opts) == 2
    assert opts[0] == opts[1] == krotov.PulseOptions(lambda_a=1.0)


def test_extract_controls(caplog):
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
        krotov.Objective(H1, X, Y),
        krotov.Objective(H1, Y, X)]
    pulse_options = {
        f: krotov.PulseOptions(lambda_a=1.0),
        g: krotov.PulseOptions(lambda_a=1.0),
    }
    controls, maps, opts, objs = krotov._extract_controls(
        objectives, pulse_options)
    assert len(controls) == 2
    assert f in controls
    assert g in controls
    assert len(maps) == len(controls)
    assert len(opts) == len(controls)
    assert len(objs) == 2

    # check error for missing PulseOptions
    pulse_options = {
        f: krotov.PulseOptions(lambda_a=1.0)}
    with pytest.raises(ValueError) as exc_info:
        krotov._extract_controls(objectives, pulse_options)
    assert 'does not have any associated pulse options' in str(exc_info.value)

    # check warning message for extra PulseOptions
    pulse_options = {
        f: krotov.PulseOptions(lambda_a=1.0),
        g: krotov.PulseOptions(lambda_a=1.0),
        h: krotov.PulseOptions(lambda_a=1.0),
    }
    with caplog.at_level(logging.WARNING):
        krotov._extract_controls(objectives, pulse_options)
    assert 'options for controls that are not in the objectives' in caplog.text
    assert len(controls) == 2

    # check same control occuring in multiple Hamiltonians
    objectives = [
        krotov.Objective(H1, X, Y),
        krotov.Objective(H2, Y, X),
        krotov.Objective(H3, Y, X)]
    pulse_options = {
        f: krotov.PulseOptions(lambda_a=1.0),
        g: krotov.PulseOptions(lambda_a=1.0),
        h: krotov.PulseOptions(lambda_a=1.0),
        d: krotov.PulseOptions(lambda_a=1.0),
    }
    controls, maps, opts, objs = krotov._extract_controls(
        objectives, pulse_options)
    assert len(controls) == 4
    for c in (f, g, h, d):
        assert c in controls
