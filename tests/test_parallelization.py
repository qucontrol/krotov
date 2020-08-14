"""Tests of krotov.parallelization."""
import io
from functools import partial

import numpy as np
import pytest
import qutip
import scipy

import krotov
import transmon_xgate_system_mod


def _func_global(x, power, *, use_sqrt=False):
    """Test function to map."""
    y = x ** power
    if use_sqrt:
        y = np.sqrt(y)
    return y


def test_parallel_map():
    """Test that parallel_map and serial_map give the same result."""

    def func(x, power, *, use_sqrt=False):
        """Test function to map."""
        # This can't be pickled so we can only map it with loky
        y = x ** power
        if use_sqrt:
            y = np.sqrt(y)
        return y

    xs = [0, 1, 2, 4]
    expected = np.sqrt(np.array(xs) ** 2)

    ys_serial = krotov.parallelization.serial_map(
        func, xs, task_args=(2,), task_kwargs={'use_sqrt': True}
    )
    assert np.max(np.abs(ys_serial - expected)) < 1e-15

    krotov.parallelization.USE_LOKY = True
    ys_parallel = krotov.parallelization.parallel_map(
        func, xs, task_args=(2,), task_kwargs={'use_sqrt': True}
    )
    assert np.max(np.abs(ys_parallel - expected)) < 1e-15

    krotov.parallelization.USE_LOKY = False
    ys_parallel = krotov.parallelization.parallel_map(
        _func_global, xs, task_args=(2,), task_kwargs={'use_sqrt': True}
    )
    assert np.max(np.abs(ys_parallel - expected)) < 1e-15

    ys_parallel = krotov.parallelization.parallel_map(
        _func_global,
        xs,
        task_args=(2,),
        task_kwargs={'use_sqrt': True},
        progress_bar=True,
    )
    assert np.max(np.abs(ys_parallel - expected)) < 1e-15


@pytest.fixture
def transmon_xgate_system():
    """Optimization system for a single-qubit transmont gate.

    See example "Optimization of an X-Gate for a Transmon Qubit".
    """

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

    tlist = np.linspace(0, 10, 100)

    H = transmon_hamiltonian()

    pulse_options = {H[1][1]: dict(lambda_a=1, update_shape=S)}

    psi0, psi1 = logical_basis(H)

    objectives = krotov.gate_objectives(
        basis_states=[psi0, psi1], gate=qutip.operators.sigmax(), H=H
    )

    return objectives, pulse_options, tlist


def test_parallel_map_fw_prop_step_loky(transmon_xgate_system):
    """Test optimization with parallel_map parameter, using loky."""
    objectives, pulse_options, tlist = transmon_xgate_system

    krotov.parallelization.USE_LOKY = True
    log = io.StringIO()
    opt_result_loky = krotov.optimize_pulses(
        objectives,
        pulse_options,
        tlist,
        propagator=krotov.propagators.expm,
        chi_constructor=krotov.functionals.chis_re,
        info_hook=partial(krotov.info_hooks.print_debug_information, out=log),
        iter_stop=1,
        skip_initial_forward_propagation=True,
        parallel_map=(
            krotov.parallelization.parallel_map,
            krotov.parallelization.parallel_map,
            krotov.parallelization.parallel_map_fw_prop_step,
        ),
    )
    tau1_loky = abs(opt_result_loky.tau_vals[0][0])
    tau2_loky = abs(opt_result_loky.tau_vals[0][1])
    assert abs(tau1_loky - 0.9693) < 1e-3
    assert abs(tau2_loky - 0.7743) < 1e-3

    krotov.parallelization.USE_LOKY = False
    # in order to be pickleable, all functions must be defined in a module
    # (transmon_xgate_system_mod)
    H = transmon_xgate_system_mod.transmon_hamiltonian()
    pulse_options = {
        H[1][1]: dict(lambda_a=1, update_shape=transmon_xgate_system_mod.S)
    }
    psi0, psi1 = transmon_xgate_system_mod.logical_basis(H)
    objectives = krotov.gate_objectives(
        basis_states=[psi0, psi1], gate=qutip.operators.sigmax(), H=H
    )
    try:
        opt_result = krotov.optimize_pulses(
            objectives,
            pulse_options,
            tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            iter_stop=1,
            skip_initial_forward_propagation=True,
            parallel_map=(
                krotov.parallelization.parallel_map,
                krotov.parallelization.parallel_map,
                krotov.parallelization.parallel_map_fw_prop_step,
            ),
        )
        tau1 = abs(opt_result.tau_vals[0][0])
        tau2 = abs(opt_result.tau_vals[0][1])
        assert abs(tau1_loky - tau1) < 1e-12
        assert abs(tau2_loky - tau2) < 1e-12
    except RuntimeError:
        # parallelization without LOKY doesn't work on all platforms
        pass
