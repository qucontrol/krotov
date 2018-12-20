"""High-level tests for `krotov` package."""

from pkg_resources import parse_version

import pytest
import qutip
import numpy as np

import krotov


def test_valid_version():
    """Check that the package defines a valid __version__"""
    assert parse_version(krotov.__version__) >= parse_version("0.0.1")


def test_complex_control_rejection():
    """Test that complex controls are rejected"""
    H0 = qutip.Qobj(0.5*np.diag([-1, 1]))
    H1 = qutip.Qobj(np.mat([[1, 2], [3, 4]]))

    psi0 = qutip.Qobj(np.array([1, 0]))
    psi1 = qutip.Qobj(np.array([0, 1]))

    def eps0(t, args):
        return 0.2 * np.exp(1j * t)

    def S(t):
        """Shape function for the field update"""
        return krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq')

    H = [H0, [H1, eps0]]

    objectives = [
        krotov.Objective(initial_state=psi0, target=psi1, H=H)
    ]

    pulse_options = {
        H[1][1]: krotov.PulseOptions(lambda_a=5, shape=S)
    }

    tlist = np.linspace(0, 5, 500)

    with pytest.raises(ValueError) as exc_info:
        krotov.optimize_pulses(
            objectives, pulse_options, tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re, iter_stop=0)
    assert 'All controls must be real-valued' in str(exc_info.value)

    def S2(t):
        """Shape function for the field update"""
        return 2.0 * krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq')


def test_reject_invalid_shapes():
    """Test that invalid control shapes are rejected"""
    H0 = qutip.Qobj(0.5*np.diag([-1, 1]))
    H1 = qutip.Qobj(np.mat([[1, 2], [3, 4]]))

    psi0 = qutip.Qobj(np.array([1, 0]))
    psi1 = qutip.Qobj(np.array([0, 1]))

    def eps0(t, args):
        return 0.2

    H = [H0, [H1, eps0]]

    objectives = [
        krotov.Objective(initial_state=psi0, target=psi1, H=H)]

    tlist = np.linspace(0, 5, 500)

    def S_complex(t):
        """Shape function for the field update"""
        return 1j * krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq')

    def S_negative(t):
        """Shape function for the field update"""
        return -1 * krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq')

    def S_large(t):
        """Shape function for the field update"""
        return 2 * krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq')

    with pytest.raises(ValueError) as exc_info:
        pulse_options = {
            H[1][1]: krotov.PulseOptions(lambda_a=5, shape=S_complex)}
        krotov.optimize_pulses(
            objectives, pulse_options, tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re, iter_stop=0)
    assert 'must be real-valued' in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        pulse_options = {
            H[1][1]: krotov.PulseOptions(lambda_a=5, shape=S_negative)}
        krotov.optimize_pulses(
            objectives, pulse_options, tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re, iter_stop=0)
    assert 'must have values in the range [0, 1]' in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        pulse_options = {
            H[1][1]: krotov.PulseOptions(lambda_a=5, shape=S_large)}
        krotov.optimize_pulses(
            objectives, pulse_options, tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re, iter_stop=0)
    assert 'must have values in the range [0, 1]' in str(exc_info.value)
