"""Tests of krotov.functionals"""
import copy
from itertools import product

import numpy as np
import pytest
import qutip
from qutip import ket

import krotov


try:
    import qutip.qip.gates as qutip_gates
except ImportError:
    import qutip.qip.operations.gates as qutip_gates


@pytest.fixture
def canonical_basis():
    return [ket('00'), ket('01'), ket('10'), ket('11')]


@pytest.fixture
def sqrt_SWAP_basis(canonical_basis):
    return krotov.functionals.mapped_basis(
        qutip_gates.sqrtswap(), canonical_basis
    )


@pytest.fixture
def cphase_objectives(canonical_basis):
    H = qutip.Qobj()  # dummy Hamiltonian (won't be used)
    return krotov.objectives.gate_objectives(
        canonical_basis, gate=qutip_gates.cphase(np.pi), H=H
    )


@pytest.fixture
def cphase_lv_full_objectives(canonical_basis):
    L = qutip.Qobj()  # dummy Liouvillian (won't be used)
    return krotov.objectives.gate_objectives(
        canonical_basis,
        gate=qutip_gates.cphase(np.pi),
        H=L,
        liouville_states_set='full',
    )


@pytest.fixture
def iswap_state_objectives(canonical_basis):
    H = qutip.Qobj()  # dummy Hamiltonian (won't be used)
    objectives = krotov.gate_objectives(
        canonical_basis, qutip_gates.sqrtiswap(), H
    )
    return objectives


@pytest.fixture
def transmon_3states_objectives():
    # see also test_objectives:test_transmon_3states_objectives
    L = qutip.Qobj()  # dummy Liouvillian (won't be used)
    n_qubit = 3
    ket00 = qutip.ket((0, 0), dim=(n_qubit, n_qubit))
    ket01 = qutip.ket((0, 1), dim=(n_qubit, n_qubit))
    ket10 = qutip.ket((1, 0), dim=(n_qubit, n_qubit))
    ket11 = qutip.ket((1, 1), dim=(n_qubit, n_qubit))
    basis = [ket00, ket01, ket10, ket11]
    weights = [20, 1, 1]
    objectives = krotov.gate_objectives(
        basis,
        qutip_gates.sqrtiswap(),
        L,
        liouville_states_set='3states',
        weights=weights,
    )
    return objectives


def test_mapped_basis_preserves_hs_structure():
    """Test that mapped_basis preserves the hilbert space structure of the
    input basis."""
    basis = [ket(nums) for nums in [(0, 0), (0, 1), (1, 0), (1, 1)]]
    states = krotov.functionals.mapped_basis(qutip_gates.cnot(), basis)
    for state in states:
        assert isinstance(state, qutip.Qobj)
        assert state.dims == basis[0].dims
        assert state.shape == basis[0].shape
        assert state.type == basis[0].type


def test_f_tau_with_weights(sqrt_SWAP_basis, cphase_objectives):
    tau_vals = [
        obj.target.overlap(psi)
        for (psi, obj) in zip(sqrt_SWAP_basis, cphase_objectives)
    ]
    assert abs(tau_vals[0] - (1 + 0j)) < 1e-14
    assert abs(tau_vals[1] - (0.5 + 0.5j)) < 1e-14
    assert abs(tau_vals[2] - (0.5 + 0.5j)) < 1e-14
    assert abs(tau_vals[3] - (-1 + 0j)) < 1e-14

    F = krotov.functionals.f_tau(sqrt_SWAP_basis, cphase_objectives)
    assert abs(F - ((1 + 1j) / 4)) < 1e-14

    objectives = copy.deepcopy(cphase_objectives)

    # objectives[0].weight = 1.0
    objectives[1].weight = 2.0
    objectives[2].weight = 0.5
    objectives[3].weight = 0
    F = krotov.functionals.f_tau(sqrt_SWAP_basis, objectives)
    assert abs(F - ((2.25 + 1.25j) / 4)) < 1e-14

    # make sure we didn't inadvertently modify the original objectives
    for obj in cphase_objectives:
        assert not hasattr(obj, 'weight')


def test_J_T_ss(sqrt_SWAP_basis, cphase_objectives):
    J = krotov.functionals.J_T_ss(sqrt_SWAP_basis, cphase_objectives)
    assert abs(J - 0.25) < 1e-14


def test_J_T_sm(sqrt_SWAP_basis, cphase_objectives):
    J = krotov.functionals.J_T_sm(sqrt_SWAP_basis, cphase_objectives)
    assert abs(J - 0.875) < 1e-14


def test_J_T_re(sqrt_SWAP_basis, cphase_objectives):
    J = krotov.functionals.J_T_re(sqrt_SWAP_basis, cphase_objectives)
    assert abs(J - 0.75) < 1e-14


def test_J_T_ss_with_weights(sqrt_SWAP_basis, cphase_objectives):
    objectives = copy.deepcopy(cphase_objectives)

    # objectives[0].weight = 1.0
    objectives[1].weight = 2.0
    objectives[2].weight = 0.5
    objectives[3].weight = 0

    J = krotov.functionals.J_T_ss(sqrt_SWAP_basis, objectives)
    assert abs(J - 1.75 / 4) < 1e-14

    for obj in cphase_objectives:
        assert not hasattr(obj, 'weight')


def test_J_T_hs_unitary(
    sqrt_SWAP_basis, cphase_objectives, cphase_lv_full_objectives
):
    """Test that for a unitary evolution, J_T_hs is equivalent to J_T_re"""
    J_hs = krotov.functionals.J_T_hs(sqrt_SWAP_basis, cphase_objectives)
    J_re = krotov.functionals.J_T_re(sqrt_SWAP_basis, cphase_objectives)
    assert abs(J_hs - J_re) < 1e-14

    # unitary density matrices should also give matching results (but not the
    # same numbers as in Hilbert space)
    J_re_1 = J_re
    fw_states_T = [
        psi * phi.dag()
        for (psi, phi) in product(sqrt_SWAP_basis, sqrt_SWAP_basis)
    ]
    J_hs = krotov.functionals.J_T_hs(fw_states_T, cphase_lv_full_objectives)
    J_re = krotov.functionals.J_T_re(fw_states_T, cphase_lv_full_objectives)
    assert abs(J_hs - J_re) < 1e-14
    assert J_re != J_re_1


def test_chi_hs_transmon(transmon_3states_objectives):
    objectives = transmon_3states_objectives
    n_qubit = objectives[0].initial_state.dims[0][0]
    ket00 = qutip.ket((0, 0), dim=(n_qubit, n_qubit))
    ket01 = qutip.ket((0, 1), dim=(n_qubit, n_qubit))
    ket10 = qutip.ket((1, 0), dim=(n_qubit, n_qubit))
    ket11 = qutip.ket((1, 1), dim=(n_qubit, n_qubit))
    ρ_mixed = 0.25 * (
        ket00 * ket00.dag()
        + ket01 * ket01.dag()
        + ket10 * ket10.dag()
        + ket11 * ket11.dag()
    )
    assert (ρ_mixed * ρ_mixed).tr() == 0.25
    assert (ρ_mixed - objectives[2].target).norm('max') < 1e-14
    fw_states_T = [ρ_mixed, ρ_mixed, ρ_mixed]
    χs = krotov.functionals.chis_hs(fw_states_T, objectives, None)
    χ1 = (1 / 6.0) * (60.0 / 22.0) * (objectives[0].target - ρ_mixed)
    χ2 = (1 / 6.0) * (3.0 / 22.0) * (objectives[1].target - ρ_mixed)
    χ3 = 0.0 * ρ_mixed
    assert (χs[0] - χ1).norm('max') < 1e-14
    assert (χs[1] - χ2).norm('max') < 1e-14
    assert (χs[2] - χ3).norm('max') < 1e-14

    # without weights
    objectives = copy.deepcopy(objectives)
    for obj in objectives:
        del obj.weight
    χs = krotov.functionals.chis_hs(fw_states_T, objectives, None)
    χ1 = (1 / 6.0) * (objectives[0].target - ρ_mixed)
    χ2 = (1 / 6.0) * (objectives[1].target - ρ_mixed)
    χ3 = 0.0 * ρ_mixed
    assert (χs[0] - χ1).norm('max') < 1e-14
    assert (χs[1] - χ2).norm('max') < 1e-14
    assert (χs[2] - χ3).norm('max') < 1e-14


def test_chi_ss(iswap_state_objectives, canonical_basis):
    χs = krotov.functionals.chis_ss(
        fw_states_T=canonical_basis,
        objectives=iswap_state_objectives,
        tau_vals=[1, 0.5 * (1 + 1j), 0.5 * (1 + 1j), 1],
    )
    χ1 = (1 / 4) * (1.0) * iswap_state_objectives[0].target
    χ2 = (1 / 4) * (0.5 * (1 + 1j)) * iswap_state_objectives[1].target
    χ3 = (1 / 4) * (0.5 * (1 + 1j)) * iswap_state_objectives[2].target
    χ4 = (1 / 4) * (1.0) * iswap_state_objectives[3].target
    assert (χs[0] - χ1).norm('max') < 1e-14
    assert (χs[1] - χ2).norm('max') < 1e-14
    assert (χs[2] - χ3).norm('max') < 1e-14
    assert (χs[3] - χ4).norm('max') < 1e-14

    # iswap_state_objectives[0].weight = 1.0
    iswap_state_objectives[1].weight = 2.0
    iswap_state_objectives[2].weight = 0.5
    iswap_state_objectives[3].weight = 0
    χs = krotov.functionals.chis_ss(
        fw_states_T=canonical_basis,
        objectives=iswap_state_objectives,
        tau_vals=[1, 0.5 * (1 + 1j), 0.5 * (1 + 1j), 1],
    )
    χ1 = (1 / 4) * 1.0 * (1.0) * iswap_state_objectives[0].target
    χ2 = (1 / 4) * 2.0 * (0.5 * (1 + 1j)) * iswap_state_objectives[1].target
    χ3 = (1 / 4) * 0.5 * (0.5 * (1 + 1j)) * iswap_state_objectives[2].target
    χ4 = (1 / 4) * 0.0 * (1.0) * iswap_state_objectives[3].target
    assert (χs[0] - χ1).norm('max') < 1e-14
    assert (χs[1] - χ2).norm('max') < 1e-14
    assert (χs[2] - χ3).norm('max') < 1e-14
    assert (χs[3] - χ4).norm('max') < 1e-14


def test_chi_sm(iswap_state_objectives, canonical_basis):
    χs = krotov.functionals.chis_sm(
        fw_states_T=canonical_basis,
        objectives=iswap_state_objectives,
        tau_vals=[1, 0.5 * (1 + 1j), 0.5 * (1 + 1j), 1],
    )
    χ1 = ((3 + 1j) / 16) * iswap_state_objectives[0].target
    χ2 = ((3 + 1j) / 16) * iswap_state_objectives[1].target
    χ3 = ((3 + 1j) / 16) * iswap_state_objectives[2].target
    χ4 = ((3 + 1j) / 16) * iswap_state_objectives[3].target
    assert (χs[0] - χ1).norm('max') < 1e-14
    assert (χs[1] - χ2).norm('max') < 1e-14
    assert (χs[2] - χ3).norm('max') < 1e-14
    assert (χs[3] - χ4).norm('max') < 1e-14

    # iswap_state_objectives[0].weight = 1.0
    iswap_state_objectives[1].weight = 2.0
    iswap_state_objectives[2].weight = 0.5
    iswap_state_objectives[3].weight = 0
    χs = krotov.functionals.chis_sm(
        fw_states_T=canonical_basis,
        objectives=iswap_state_objectives,
        tau_vals=[1, 0.5 * (1 + 1j), 0.5 * (1 + 1j), 1],
    )
    χ1 = ((2.25 + 1.25j) / 16) * 1.0 * iswap_state_objectives[0].target
    χ2 = ((2.25 + 1.25j) / 16) * 2.0 * iswap_state_objectives[1].target
    χ3 = ((2.25 + 1.25j) / 16) * 0.5 * iswap_state_objectives[2].target
    χ4 = ((2.25 + 1.25j) / 16) * 0.0 * iswap_state_objectives[3].target
    assert (χs[0] - χ1).norm('max') < 1e-14
    assert (χs[1] - χ2).norm('max') < 1e-14
    assert (χs[2] - χ3).norm('max') < 1e-14
    assert (χs[3] - χ4).norm('max') < 1e-14


def test_chi_re(iswap_state_objectives, canonical_basis):
    χs = krotov.functionals.chis_re(
        canonical_basis, iswap_state_objectives, None
    )
    χ1 = (1 / 8) * iswap_state_objectives[0].target
    χ2 = (1 / 8) * iswap_state_objectives[1].target
    χ3 = (1 / 8) * iswap_state_objectives[2].target
    χ4 = (1 / 8) * iswap_state_objectives[3].target
    assert (χs[0] - χ1).norm('max') < 1e-14
    assert (χs[1] - χ2).norm('max') < 1e-14
    assert (χs[2] - χ3).norm('max') < 1e-14
    assert (χs[3] - χ4).norm('max') < 1e-14

    # iswap_state_objectives[0].weight = 1.0
    iswap_state_objectives[1].weight = 2.0
    iswap_state_objectives[2].weight = 0.5
    iswap_state_objectives[3].weight = 0
    χs = krotov.functionals.chis_re(
        canonical_basis, iswap_state_objectives, None
    )
    χ1 = (1 / 8) * 1.0 * iswap_state_objectives[0].target
    χ2 = (1 / 8) * 2.0 * iswap_state_objectives[1].target
    χ3 = (1 / 8) * 0.5 * iswap_state_objectives[2].target
    χ4 = (1 / 8) * 0.0 * iswap_state_objectives[3].target
    assert (χs[0] - χ1).norm('max') < 1e-14
    assert (χs[1] - χ2).norm('max') < 1e-14
    assert (χs[2] - χ3).norm('max') < 1e-14
    assert (χs[3] - χ4).norm('max') < 1e-14


def test_F_avg_psi(sqrt_SWAP_basis, canonical_basis):
    F = krotov.functionals.F_avg(
        fw_states_T=sqrt_SWAP_basis,
        basis_states=canonical_basis,
        gate=qutip_gates.cphase(np.pi),
    )
    assert abs(F - 0.3) < 1e-14


def test_F_avg_rho(sqrt_SWAP_basis, canonical_basis):
    fw_states_T = [
        psi * phi.dag()
        for (psi, phi) in product(sqrt_SWAP_basis, sqrt_SWAP_basis)
    ]
    F = krotov.functionals.F_avg(
        fw_states_T=fw_states_T,
        basis_states=canonical_basis,
        gate=qutip_gates.cphase(np.pi),
    )
    assert abs(F - 0.3) < 1e-14
