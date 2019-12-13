"""High-level tests for `krotov` package."""

import io
import logging
import os
from copy import deepcopy
from shutil import copyfile

import numpy as np
import pytest
import qutip
from pkg_resources import parse_version

import krotov


def test_valid_version():
    """Check that the package defines a valid __version__"""
    assert parse_version(krotov.__version__) >= parse_version("0.1.0")


def test_complex_control_rejection():
    """Test that complex controls are rejected"""
    H0 = qutip.Qobj(0.5 * np.diag([-1, 1]))
    H1 = qutip.Qobj(np.mat([[1, 2], [3, 4]]))

    psi0 = qutip.Qobj(np.array([1, 0]))
    psi1 = qutip.Qobj(np.array([0, 1]))

    def eps0(t, args):
        return 0.2 * np.exp(1j * t)

    def S(t):
        """Shape function for the field update"""
        return krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq'
        )

    H = [H0, [H1, eps0]]

    objectives = [krotov.Objective(initial_state=psi0, target=psi1, H=H)]

    pulse_options = {H[1][1]: dict(lambda_a=5, update_shape=S)}

    tlist = np.linspace(0, 5, 500)

    with pytest.raises(ValueError) as exc_info:
        krotov.optimize_pulses(
            objectives,
            pulse_options,
            tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            iter_stop=0,
        )
    assert 'all controls must be real-valued' in str(exc_info.value)

    def S2(t):
        """Shape function for the field update"""
        return 2.0 * krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq'
        )


def test_reject_invalid_shapes():
    """Test that invalid control shapes are rejected"""
    H0 = qutip.Qobj(0.5 * np.diag([-1, 1]))
    H1 = qutip.Qobj(np.mat([[1, 2], [3, 4]]))

    psi0 = qutip.Qobj(np.array([1, 0]))
    psi1 = qutip.Qobj(np.array([0, 1]))

    def eps0(t, args):
        return 0.2

    H = [H0, [H1, eps0]]

    objectives = [krotov.Objective(initial_state=psi0, target=psi1, H=H)]

    tlist = np.linspace(0, 5, 500)

    def S_complex(t):
        """Shape function for the field update"""
        return 1j * krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq'
        )

    def S_negative(t):
        """Shape function for the field update"""
        return -1 * krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq'
        )

    def S_large(t):
        """Shape function for the field update"""
        return 2 * krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq'
        )

    with pytest.raises(ValueError) as exc_info:
        pulse_options = {H[1][1]: dict(lambda_a=5, update_shape=S_complex)}
        krotov.optimize_pulses(
            objectives,
            pulse_options,
            tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            iter_stop=0,
        )
    assert 'must be real-valued' in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        pulse_options = {H[1][1]: dict(lambda_a=5, update_shape=S_negative)}
        krotov.optimize_pulses(
            objectives,
            pulse_options,
            tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            iter_stop=0,
        )
    assert 'must have values in the range [0, 1]' in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        pulse_options = {H[1][1]: dict(lambda_a=5, update_shape=S_large)}
        krotov.optimize_pulses(
            objectives,
            pulse_options,
            tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            iter_stop=0,
        )
    assert 'must have values in the range [0, 1]' in str(exc_info.value)


@pytest.fixture
def simple_state_to_state_system():
    """System from 01_example_simple_state_to_state.ipynb"""
    omega = 1.0
    ampl0 = 0.2

    H0 = -0.5 * omega * qutip.operators.sigmaz()
    H1 = qutip.operators.sigmax()
    eps0 = lambda t, args: ampl0
    H = [H0, [H1, eps0]]

    psi0 = qutip.ket('0')
    psi1 = qutip.ket('1')

    objectives = [krotov.Objective(initial_state=psi0, target=psi1, H=H)]

    def S(t):
        """Shape function for the field update"""
        return krotov.shapes.flattop(
            t, t_start=0, t_stop=5, t_rise=0.3, t_fall=0.3, func='sinsq'
        )

    pulse_options = {H[1][1]: dict(lambda_a=5, update_shape=S)}

    tlist = np.linspace(0, 5, 500)

    return objectives, pulse_options, tlist


@pytest.mark.parametrize('iter_stop', [0, -1])
def test_zero_iterations(iter_stop, simple_state_to_state_system):
    objectives, pulse_options, tlist = simple_state_to_state_system

    with io.StringIO() as log_fh:

        result = krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            store_all_pulses=True,
            info_hook=krotov.info_hooks.print_table(
                J_T=krotov.functionals.J_T_re, out=log_fh
            ),
            iter_stop=iter_stop,
            skip_initial_forward_propagation=True,
        )

        log = log_fh.getvalue()

    assert len(log.splitlines()) == 2
    assert result.message == 'Reached 0 iterations'
    assert len(result.guess_controls) == 1  # one control field
    assert len(result.guess_controls) == len(result.optimized_controls)
    assert len(result.guess_controls[0]) == len(result.tlist)
    assert len(result.optimized_controls[0]) == len(result.tlist)
    for (c1, c2) in zip(result.guess_controls, result.optimized_controls):
        assert np.all(c1 == c2)
    for pulses_for_iteration in result.all_pulses:
        for pulse in pulses_for_iteration:
            # the pulses are defined on the *intervals* of tlist
            assert len(pulse) == len(result.tlist) - 1


def test_continue_optimization(
    simple_state_to_state_system, request, tmpdir, caplog
):
    """Big integration test for a simple optimization, with various
    uses of `continue_from`. This covers a lot of edge cases not covered by the
    example notebooks.
    """
    objectives, pulse_options, tlist = simple_state_to_state_system

    dumpfile = str(tmpdir.join("oct_result_{iter:03d}.dump"))
    logfile = str(tmpdir.join("oct.log"))
    testdir = os.path.splitext(request.module.__file__)[0]
    logfile_expected = os.path.join(testdir, 'oct.log')

    with open(logfile, 'w', encoding='utf8') as log_fh:

        # initial optimization
        with caplog.at_level(logging.WARNING):
            oct_result1 = krotov.optimize_pulses(
                objectives,
                pulse_options=pulse_options,
                tlist=tlist,
                propagator=krotov.propagators.expm,
                chi_constructor=krotov.functionals.chis_re,
                store_all_pulses=True,
                info_hook=krotov.info_hooks.print_table(
                    J_T=krotov.functionals.J_T_re, out=log_fh
                ),
                check_convergence=krotov.convergence.Or(
                    krotov.convergence.check_monotonic_error,
                    krotov.convergence.dump_result(dumpfile, every=2),
                ),
                iter_stop=3,
                skip_initial_forward_propagation=True,
                # not officially supported, but should work in this case
            )
        assert (
            "You should not use `skip_initial_forward_propagation`"
            in caplog.text
        )

        # fmt: off
        assert len(oct_result1.iters) == 4  # 0 ... 3
        assert len(oct_result1.iter_seconds) == 4
        assert len(oct_result1.info_vals) == 4
        assert len(oct_result1.all_pulses) == 4
        assert len(oct_result1.states) == 1
        assert len(oct_result1.guess_controls) == 1  # one control field
        assert len(oct_result1.guess_controls) == len(oct_result1.optimized_controls)
        assert len(oct_result1.guess_controls[0]) == len(oct_result1.tlist)
        assert len(oct_result1.optimized_controls[0]) == len(oct_result1.tlist)
        for pulses_for_iteration in oct_result1.all_pulses:
            for pulse in pulses_for_iteration:
                # the pulses are defined on the *intervals* of tlist
                assert len(pulse) == len(oct_result1.tlist) - 1
        assert "3 iterations" in oct_result1.message
        # fmt: on

        # repeating the same optimization only propagates the guess pulse
        # (we'll check this later while verifying the output of the log file)
        krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            store_all_pulses=True,
            info_hook=krotov.info_hooks.print_table(
                J_T=krotov.functionals.J_T_re, out=log_fh
            ),
            check_convergence=krotov.convergence.Or(
                krotov.convergence.check_monotonic_error,
                krotov.convergence.dump_result(dumpfile, every=2),
            ),
            continue_from=oct_result1,
            iter_stop=3,
        )

        # another 2 iterations
        oct_result2 = krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            store_all_pulses=True,
            info_hook=krotov.info_hooks.print_table(
                J_T=krotov.functionals.J_T_re, out=log_fh
            ),
            check_convergence=krotov.convergence.Or(
                krotov.convergence.check_monotonic_error,
                krotov.convergence.dump_result(dumpfile, every=2),
            ),
            continue_from=oct_result1,
            iter_stop=5,
        )

        assert len(oct_result2.iters) == 6  # 0 ... 5
        assert len(oct_result2.iter_seconds) == 6
        assert len(oct_result2.info_vals) == 6
        assert len(oct_result2.all_pulses) == 6
        assert len(oct_result2.states) == 1
        assert "5 iterations" in oct_result2.message

        # and 2 more (skipping initial propagation)
        oct_result3 = krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            store_all_pulses=True,
            info_hook=krotov.info_hooks.print_table(
                J_T=krotov.functionals.J_T_re, out=log_fh
            ),
            check_convergence=krotov.convergence.Or(
                krotov.convergence.check_monotonic_error,
                krotov.convergence.dump_result(dumpfile, every=2),
            ),
            continue_from=oct_result2,
            iter_stop=7,
            skip_initial_forward_propagation=True,
        )

        assert len(oct_result3.iters) == 8  # 0 ... 7
        assert len(oct_result3.iter_seconds) == 8
        assert len(oct_result3.info_vals) == 8
        assert len(oct_result3.all_pulses) == 8
        assert len(oct_result3.states) == 1
        assert "7 iterations" in oct_result3.message

        # no-op: already anough iterations
        oct_result4 = krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            store_all_pulses=True,
            info_hook=krotov.info_hooks.print_table(
                J_T=krotov.functionals.J_T_re, out=log_fh
            ),
            check_convergence=krotov.convergence.Or(
                krotov.convergence.check_monotonic_error,
                krotov.convergence.dump_result(dumpfile, every=2),
            ),
            continue_from=oct_result3,
            iter_stop=5,  # < 7  in oct_result3
            skip_initial_forward_propagation=True,
        )

        assert len(oct_result4.iters) == 8  # 0 ... 7
        assert len(oct_result4.iter_seconds) == 8
        assert len(oct_result4.info_vals) == 8
        assert len(oct_result4.all_pulses) == 8
        assert len(oct_result4.states) == 1
        assert "7 iterations" in oct_result4.message

        assert (
            oct_result4.start_local_time_str
            == oct_result3.start_local_time_str
        )
        assert oct_result4.iters == oct_result3.iters
        assert oct_result4.message == oct_result3.message

    # fmt: off
    # check the combined log file
    if not os.path.isfile(logfile_expected):
        copyfile(logfile, logfile_expected)
    with open(logfile, encoding='utf8') as fh1, open(logfile_expected, encoding='utf8') as fh2:
        for line1 in fh1:
            line2 = next(fh2)
            assert line1[:63] == line2[:63]
    # fmt: on

    # continue from an incomplete dump file

    with caplog.at_level(logging.WARNING):
        result = krotov.result.Result.load(
            str(tmpdir.join("oct_result_004.dump"))
        )
    assert 'Result.objectives contains control placeholders' in caplog.text
    with pytest.raises(ValueError) as exc_info:
        oct_result5 = krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            store_all_pulses=True,
            check_convergence=krotov.convergence.Or(
                krotov.convergence.check_monotonic_error,
                krotov.convergence.dump_result(dumpfile, every=2),
            ),
            continue_from=result,
            iter_stop=7,
            skip_initial_forward_propagation=True,
        )
    assert "objectives must remain unchanged" in str(exc_info.value)
    result = krotov.result.Result.load(
        str(tmpdir.join("oct_result_004.dump")), objectives=objectives
    )
    assert result.iters[-1] == 4
    oct_result5 = krotov.optimize_pulses(
        objectives,
        pulse_options=pulse_options,
        tlist=tlist,
        propagator=krotov.propagators.expm,
        chi_constructor=krotov.functionals.chis_re,
        store_all_pulses=True,
        info_hook=krotov.functionals.J_T_re,
        check_convergence=krotov.convergence.Or(
            krotov.convergence.check_monotonic_error,
            krotov.convergence.dump_result(dumpfile, every=2),
        ),
        continue_from=result,
        iter_stop=7,
        skip_initial_forward_propagation=True,
    )
    assert oct_result5.iters == oct_result3.iters
    assert len(oct_result5.iter_seconds) == 8
    assert len(oct_result5.info_vals) == 8
    assert len(oct_result5.all_pulses) == 8
    assert "7 iterations" in oct_result5.message
    Δ = np.max(
        np.abs(
            oct_result5.optimized_controls[-1]
            - oct_result3.optimized_controls[-1]
        )
    )
    assert Δ < 1e-10

    # broken continuation: different number of objectives
    result_with_extra_objective = deepcopy(result)
    result_with_extra_objective.objectives.append(
        deepcopy(result.objectives[0])
    )
    with pytest.raises(ValueError) as exc_info:
        krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            continue_from=result_with_extra_objective,
        )
    assert "number of objectives must be the same" in str(exc_info.value)

    # broken continuation: different value for store_all_pulses
    with pytest.raises(ValueError) as exc_info:
        krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            continue_from=result,
            store_all_pulses=False,
        )
    assert "store_all_pulses parameter cannot be changed" in str(
        exc_info.value
    )

    # broken continuation: different time units
    result_with_scaled_tlist = deepcopy(result)
    result_with_scaled_tlist.objectives = result.objectives
    result_with_scaled_tlist.tlist *= 2
    with pytest.raises(ValueError) as exc_info:
        krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            continue_from=result_with_scaled_tlist,
            store_all_pulses=True,
        )
    assert "same time grid" in str(exc_info.value)

    # broken continuation: changed nt
    result_with_changed_nt = deepcopy(result)
    result_with_changed_nt.objectives = result.objectives
    result_with_changed_nt.tlist = np.linspace(0, 5, 1000)
    with pytest.raises(ValueError) as exc_info:
        krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            continue_from=result_with_changed_nt,
            store_all_pulses=True,
        )
    assert "same time grid" in str(exc_info.value)

    # incongruent controls
    result_with_incongruent_pulse = deepcopy(result)
    result_with_incongruent_pulse.objectives = result.objectives
    result_with_incongruent_pulse.optimized_controls[0] = np.stack(
        [result.optimized_controls[0], result.optimized_controls[0]]
    ).flatten()
    with pytest.raises(ValueError) as exc_info:
        krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            continue_from=result_with_incongruent_pulse,
            store_all_pulses=True,
        )
    assert "optimized_controls and tlist are incongruent" in str(
        exc_info.value
    )

    # passing complete garbage to `continue_from`
    with pytest.raises(ValueError) as exc_info:
        krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            continue_from=result.objectives[0],
            store_all_pulses=True,
        )
    assert "only possible from a Result object" in str(exc_info.value)
