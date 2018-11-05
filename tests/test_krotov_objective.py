"""Tests for KrotovObjective in isolation"""
import numpy as np
import scipy
import qutip

from krotov import KrotovObjective

import pytest


@pytest.fixture
def transmon_ham_and_states(
        Ec=0.386, EjEc=45, nstates=8, ng=0.0, T=10.0):
    """Transmon Hamiltonian"""
    Ej = EjEc * Ec
    n = np.arange(-nstates, nstates+1)
    up = np.diag(np.ones(2*nstates), k=-1)
    do = up.T
    H0 = qutip.Qobj(np.diag(4*Ec*(n - ng)**2) - Ej*(up+do)/2.0)
    H1 = qutip.Qobj(-2*np.diag(n))

    eigenvals, eigenvecs = scipy.linalg.eig(H0.full())
    ndx = np.argsort(eigenvals.real)
    E = eigenvals[ndx].real
    V = eigenvecs[:, ndx]
    w01 = E[1]-E[0]  # Transition energy between states

    psi0 = qutip.Qobj(V[:, 0])
    psi1 = qutip.Qobj(V[:, 1])

    profile = lambda t: np.exp(-40.0*(t/T - 0.5)**2)
    eps0 = lambda t, args: 0.5 * profile(t) * np.cos(8*np.pi*w01*t)
    return ([H0, [H1, eps0]], psi0, psi1)


def test_krotov_objective_initialization(transmon_ham_and_states):
    """Test basic instantiation of a KrotovObjective with qutip objects"""
    H, psi0, psi1 = transmon_ham_and_states
    target = KrotovObjective(H, psi0, psi1)
    assert target.H == H
    assert target.initial_state == psi0
    assert target.target_state == psi1
    assert target == KrotovObjective(
        H=H, initial_state=psi0, target_state=psi1)
