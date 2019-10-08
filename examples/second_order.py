import functools

import numpy as np
import qutip
import weylchamber as wc  # pip install weylchamber
from weylchamber.coordinates import from_magic

import krotov


w1 = 1.1  # qubit 1 level splitting
w2 = 2.1  # qubit 2 level splitting
J = 0.2  # effective qubit coupling
u0 = 0.3  # initial driving strength
la = 1.1  # relative pulse coupling strength of second qubit
T = 25.0  # final time
nt = 250  # number of time steps

tlist = np.linspace(0, T, nt)


def eps0(t, args):
    return u0 * krotov.shapes.flattop(
        t, t_start=0, t_stop=T, t_rise=(T / 20), t_fall=(T / 20), func='sinsq'
    )


def hamiltonian(w1=w1, w2=w2, J=J, la=la, u0=u0):
    """Two qubit Hamiltonian

    Args:
        w1 (float): energy separation of the first qubit levels
        w2 (float): energy separation of the second qubit levels
        J (float): effective coupling between both qubits
        la (float): factor that pulse coupling strength differs for second
            qubit
        u0 (float): constant amplitude of the driving field
    """
    # local qubit Hamiltonians
    Hq1 = 0.5 * w1 * np.diag([-1, 1])
    Hq2 = 0.5 * w2 * np.diag([-1, 1])

    # lift Hamiltonians to joint system operators
    H0 = np.kron(Hq1, np.identity(2)) + np.kron(np.identity(2), Hq2)

    # define the interaction Hamiltonian
    sig_x = np.array([[0, 1], [1, 0]])
    sig_y = np.array([[0, -1j], [1j, 0]])
    Hint = 2 * J * (np.kron(sig_x, sig_x) + np.kron(sig_y, sig_y))
    H0 = H0 + Hint

    # define the drive Hamiltonian
    H1 = np.kron(np.array([[0, 1], [1, 0]]), np.identity(2)) + la * np.kron(
        np.identity(2), np.array([[0, 1], [1, 0]])
    )

    # convert Hamiltonians to QuTiP objects
    H0 = qutip.Qobj(H0)
    H1 = qutip.Qobj(H1)

    return [H0, [H1, eps0]]


H = hamiltonian(w1=w1, w2=w2, J=J, la=la, u0=u0)

psi_00 = qutip.Qobj(np.kron(np.array([1, 0]), np.array([1, 0])))
psi_01 = qutip.Qobj(np.kron(np.array([1, 0]), np.array([0, 1])))
psi_10 = qutip.Qobj(np.kron(np.array([0, 1]), np.array([1, 0])))
psi_11 = qutip.Qobj(np.kron(np.array([0, 1]), np.array([0, 1])))

objectives = krotov.gate_objectives(
    basis_states=[psi_00, psi_01, psi_10, psi_11], gate="PE", H=H
)

chi_constructor = wc.perfect_entanglers.make_PE_krotov_chi_constructor(
    [psi_00, psi_01, psi_10, psi_11]
)


def S(t):
    """Shape function for the field update"""
    return krotov.shapes.flattop(
        t, t_start=0, t_stop=T, t_rise=T / 20, t_fall=T / 20, func='sinsq'
    )


pulse_options = {H[1][1]: dict(lambda_a=1.0e2, update_shape=S)}


def calculate_PE_val(**args):
    basis = [objectives[i].initial_state for i in [0, 1, 2, 3]]
    states = [args['fw_states_T'][i] for i in [0, 1, 2, 3]]
    U = wc.gates.gate(basis, states)
    c1, c2, c3 = wc.coordinates.c1c2c3(from_magic(U))
    g1, g2, g3 = wc.local_invariants.g1g2g3_from_c1c2c3(c1, c2, c3)
    conc = wc.perfect_entanglers.concurrence(c1, c2, c3)
    F_PE = wc.perfect_entanglers.F_PE(g1, g2, g3)
    print("    F_PE: %.2f\n    gate conc.: %.2f" % (F_PE, conc))
    return F_PE


# we don't really need this example to do the actual optimization, so we'll
# just stop after the first iteration, to make this run through quickly
krotov.optimize_pulses = functools.partial(krotov.optimize_pulses, iter_stop=1)

# fmt: off
###############################################################################
class sigma(krotov.second_order.Sigma):
    def __init__(self, A, epsA=0):
        self.A = A
        self.epsA = epsA

    def __call__(self, t):
        return -max(self.epsA, 2 * self.A + self.epsA)

    def refresh(
        self, forward_states, forward_states0,
        chi_states, chi_norms, optimized_pulses,
        guess_pulses, objectives, result,
    ):
        try:
            # info_vals contains values of PE functional
            Delta_J_T = (
                result.info_vals[-1][0] - result.info_vals[-2][0]
            )
        except IndexError:  # first iteration
            Delta_J_T = 0
        self.A = krotov.second_order.numerical_estimate_A(
            forward_states, forward_states0, chi_states,
            chi_norms, Delta_J_T
        )


oct_result = krotov.optimize_pulses(
    objectives,
    pulse_options=pulse_options,
    tlist=tlist,
    propagator=krotov.propagators.expm,
    chi_constructor=chi_constructor,  # from weylchamber package
    info_hook=calculate_PE_val,
    sigma=sigma(A=0.0),
)
###############################################################################
# fmt: on
