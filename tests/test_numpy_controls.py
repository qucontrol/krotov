"""Test using numpy controls."""
import numpy as np
import pytest
import qutip

import krotov


@pytest.fixture
def numpy_control_system():
    """Optimization system for State-to-State Transfer

    This is the same as in 01_example_simple_state_to_state.ipynb, but using a
    numpy array as the control.
    """
    tlist = np.linspace(0, 5, 500)

    def hamiltonian(omega=1.0, ampl0=0.2):
        """Two-level-system Hamiltonian

        Args:
            omega (float): energy separation of the qubit levels
            ampl0 (float): constant amplitude of the driving field
        """
        H0 = -0.5 * omega * qutip.operators.sigmaz()
        H1 = qutip.operators.sigmax()

        def guess_control(t, args):
            return ampl0 * krotov.shapes.flattop(
                t, t_start=0, t_stop=5, t_rise=0.3, func="blackman"
            )

        guess_array = np.array([guess_control(t, []) for t in tlist])
        # return [H0, [H1, guess_control]]
        return [H0, [H1, guess_array]]

    H = hamiltonian()

    objectives = [
        krotov.Objective(
            initial_state=qutip.ket("0"),
            target=qutip.ket("1"),
            H=H,
        )
    ]

    def S(t):
        """Shape function for the field update"""
        return krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='blackman'
        )

    # pulse_options = {H[1][1]: dict(lambda_a=5, update_shape=S)}
    pulse_options = {id(H[1][1]): dict(lambda_a=5, update_shape=S)}

    return tlist, objectives, pulse_options


def test_numpy_controls(numpy_control_system):
    """Test optimization with numpy array controls.

    Test the resolution of https://github.com/qucontrol/krotov/issues/79
    """
    tlist, objectives, pulse_options = numpy_control_system
    opt_result = krotov.optimize_pulses(
        objectives,
        pulse_options=pulse_options,
        tlist=tlist,
        propagator=krotov.propagators.expm,
        chi_constructor=krotov.functionals.chis_ss,
        check_convergence=krotov.convergence.Or(
            krotov.convergence.value_below('1e-3', name='J_T'),
            krotov.convergence.check_monotonic_error,
        ),
        iter_stop=0,
        skip_initial_forward_propagation=True,
    )
    assert opt_result.message == 'Reached 0 iterations'
