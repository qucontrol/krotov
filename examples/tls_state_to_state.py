#!/usr/bin/env python
"""Example script for the optimization of a simple state-to-state
transition in a two-level system"""
import krotov
import qutip
import numpy as np


# First, we define the physical system (a simple TLS)

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

    return [H0, [H1, guess_control]]


H = hamiltonian()
tlist = np.linspace(0, 5, 500)


# Second, we define the control objective: a state-to-state
# transition from the |0⟩ eigenstate to the |1⟩ eigenstate

objectives = [
    krotov.Objective(
        initial_state=qutip.ket("0"), target=qutip.ket("1"), H=H
    )
]


# The magnitude of the pulse updates at each point in time are
# determined by the Krotov step size lambda_a and the
# time-dependent update shape (in [0, 1])
def S(t):
    """Shape function for the field update"""
    return krotov.shapes.flattop(
        t, t_start=0, t_stop=5, t_rise=0.3, func="blackman"
    )


# set required parameters for H[1][1] (the guess_control)
pulse_options = {H[1][1]: dict(lambda_a=5, update_shape=S)}


# Before performing the optimization, it is usually a good idea
# to observe the system dynamics under the guess pulse. The
# mesolve method of the objective delegates to QuTiP's mesolve,
# and can calculate the expectation values of the projectors
# onto the |0⟩ and |1⟩ states, i.e., the population.

proj0, proj1 = (qutip.ket2dm(qutip.ket(l)) for l in ("0", "1"))
e_ops = [proj0, proj1]
guess_dynamics = objectives[0].mesolve(tlist, e_ops=e_ops)

# the resulting expectations values are in guess_dynamics.expect.
# The final-time populations are:

print(
    "guess final time population in |0⟩, |1⟩: %.3f, %.3f\n"
    % tuple([guess_dynamics.expect[l][-1] for l in (0, 1)])
)


# Now, we perform the actual optimization

opt_result = krotov.optimize_pulses(
    objectives,
    pulse_options=pulse_options,
    tlist=tlist,
    propagator=krotov.propagators.expm,
    chi_constructor=krotov.functionals.chis_ss,
    info_hook=krotov.info_hooks.print_table(
        J_T=krotov.functionals.J_T_ss
    ),
    check_convergence=krotov.convergence.Or(
        krotov.convergence.value_below('1e-3', name='J_T'),
        krotov.convergence.check_monotonic_error,
    ),
    store_all_pulses=True,
)

print("\n", opt_result, sep='')


# We can observe the population dynamics under the optimized
# control

opt_dynamics = opt_result.optimized_objectives[0].mesolve(
    tlist, e_ops=[proj0, proj1]
)

print(
    "\noptimized final time population in |0⟩, |1⟩: %.3f, %.3f"
    % (opt_dynamics.expect[0][-1], opt_dynamics.expect[1][-1])
)
#################################################################
import os

outfile = os.path.splitext(os.path.abspath(__file__))[0] + '.dump'
opt_result.dump(outfile)
# vim: set textwidth=65: (SciPost formatting)
