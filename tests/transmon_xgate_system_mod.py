"""The transmon_xgate_system fixture for test_parallelization.py in module
form.

This needs to be in a module so that all the functions are pickleable

"""
import numpy as np
import qutip
import scipy

import krotov


def eps0(t, args):
    T = 10
    return 4 * np.exp(-40.0 * (t / T - 0.5) ** 2)


def transmon_hamiltonian(Ec=0.386, EjEc=45, nstates=2, ng=0.0, T=10.0):
    Ej = EjEc * Ec
    n = np.arange(-nstates, nstates + 1)
    up = np.diag(np.ones(2 * nstates), k=-1)
    do = up.T
    H0 = qutip.Qobj(np.diag(4 * Ec * (n - ng) ** 2) - Ej * (up + do) / 2.0)
    H1 = qutip.Qobj(-2 * np.diag(n))

    return [H0, [H1, eps0]]


def logical_basis(H):
    H0 = H[0]
    eigenvals, eigenvecs = scipy.linalg.eig(H0.full())
    ndx = np.argsort(eigenvals.real)
    V = eigenvecs[:, ndx]
    psi0 = qutip.Qobj(V[:, 0])
    psi1 = qutip.Qobj(V[:, 1])
    return psi0, psi1


def S(t):
    return krotov.shapes.flattop(
        t, t_start=0.0, t_stop=10.0, t_rise=0.5, func='sinsq'
    )
