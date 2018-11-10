"""Tests for krotov.Objective in isolation"""
import copy

import numpy as np
import scipy
import qutip

import krotov

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
    """Test basic instantiation of a krotov.Objective with qutip objects"""
    H, psi0, psi1 = transmon_ham_and_states
    target = krotov.Objective(H, psi0, psi1)
    assert target.H == H
    assert target.initial_state == psi0
    assert target.target_state == psi1
    assert target == krotov.Objective(
        H=H, initial_state=psi0, target_state=psi1)


def test_objective_copy(transmon_ham_and_states):
    """Test that copy.copy(objective) produces the expected equalities by value
    and by reference"""
    H, psi0, psi1 = transmon_ham_and_states
    c1 = H[1].copy()  # we just need something structurally sound ...
    c2 = H[1].copy()  # ... It doesn't need to make sense physically
    assert c1 == c2      # equal by value
    assert c1 is not c2  # not equal by reference

    target1 = krotov.Objective(H, psi0, psi1, c_ops=[c1, c2])
    target2 = copy.copy(target1)
    assert target1 == target2
    assert target1 is not target2
    assert target1.H == target2.H
    assert target1.H is not target2.H
    assert target1.H[0] is target2.H[0]
    assert target1.H[1] is not target2.H[1]
    assert target1.H[1][0] is target2.H[1][0]
    assert target1.H[1][1] is target2.H[1][1]
    assert target1.c_ops[0] == target2.c_ops[0]
    assert target1.c_ops[0] is not target2.c_ops[0]
    assert target1.c_ops[0][0] is target2.c_ops[0][0]
    assert target1.c_ops[0][1] is target2.c_ops[0][1]
