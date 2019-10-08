import qutip
import krotov


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
            t, t_start=0, t_stop=5, t_rise=0.3, func="sinsq"
        )

    return [H0, [H1, guess_control]]


H = hamiltonian()

###############################################################################
objectives = krotov.gate_objectives(
    basis_states=[qutip.ket('0'), qutip.ket('1')],
    gate=qutip.operators.sigmax(),
    H=H,
)
###############################################################################

krotov.objectives.Objective.str_use_unicode = False
print(objectives)
