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


def L(omega, delta):
    """Two-level-system Liouvillian

    Args:
        omega (float): energy separation of the qubit levels
        delta (float): factor in front of the default amplitude (0.2)
    """
    A = qutip.operators.sigmam()
    H = hamiltonian(omega, ampl0=delta*0.2)
    return krotov.objectives.liouvillian(H, c_ops=[A])


omega_vals = [0.9, 0.95, 1.05, 1.1]
delta_vals = [0.9, 0.95, 1.05, 1.1]


###############################################################################
from math import pi  # standard library
objectives = krotov.gate_objectives(
    basis_states=[
        qutip.ket(l) for l in ['00', '01', '10', '11']
    ],
    gate=qutip.gates.cphase(pi),
    H=L(omega=1, delta=0),
    liouville_states_set='3states',
    weights=[0, 1, 1]
)
###############################################################################
import itertools  # standard library
ensemble_liouvillians = [
    L(omega, delta)
    for (omega, delta)
    in itertools.product(omega_vals, delta_vals)
]
objectives = krotov.objectives.ensemble_objectives(
    objectives, ensemble_liouvillians
)
###############################################################################

krotov.objectives.Objective.str_use_unicode = False
print(objectives)
