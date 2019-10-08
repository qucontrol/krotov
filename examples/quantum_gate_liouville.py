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


A = qutip.operators.sigmam()
H = hamiltonian()
L = krotov.objectives.liouvillian(H, c_ops=[A])

###############################################################################
objectives = krotov.gate_objectives(
    basis_states=[qutip.ket(l) for l in ['00', '01', '10', '11']],
    gate=qutip.gates.sqrtiswap(),
    H=L,  # Liouvillian super-operator (qutip.Qobj instance)
    liouville_states_set='3states',
    weights=[20, 1, 1],
)
###############################################################################

krotov.objectives.Objective.str_use_unicode = False
print(objectives)
