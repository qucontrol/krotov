"""Test for modify_params_after_iter/info_hook"""
import contextlib
import io
import os
from pathlib import Path

import numpy as np
import pytest

import krotov
from krotov.info_hooks import print_table
from test_krotov import simple_state_to_state_system
from test_objectives import transmon_ham_and_states


def test_infohook_chaining(transmon_ham_and_states):
    """Test that transmon_ham_and_states and info_hooks get chained together
    correctly. This tests a whole bunch of implementation details:

    - return values from multiple info_hooks combine in tuple
    - return value None (from modify_params_after_iter) gets ignored
    - shared_data gets passed along through multiple hooks
    - shared_data is cleared in each iteration
    """
    H, psi0, psi1 = transmon_ham_and_states
    obj = krotov.Objective(initial_state=psi0, target=psi1, H=H)
    tlist = np.array([0, 0.01, 0.02])
    pulse_options = {H[1][1]: dict(lambda_a=1, update_shape=1)}
    stdout = io.StringIO()

    def adjust_lambda_a(**args):
        λₐ = args['lambda_vals'][0]
        args['lambda_vals'][0] *= 0.5
        if 'messages' not in args['shared_data']:
            args['shared_data']['messages'] = []
        args['shared_data']['messages'].append(
            'λₐ: %s → %s' % (λₐ, args['lambda_vals'][0])
        )

    def print_fidelity(**args):
        F_re = np.average(np.array(args['tau_vals']).real)
        print("Iteration %d: \tF = %f" % (args['iteration'], F_re))
        return F_re

    def print_messages(**args):
        if 'messages' in args['shared_data']:
            message = "; ".join(
                [msg for msg in args['shared_data']['messages']]
            )
            print("\tmsg: " + message)
            return message

    with contextlib.redirect_stdout(stdout):
        oct_result = krotov.optimize_pulses(
            [obj],
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            info_hook=krotov.info_hooks.chain(print_fidelity, print_messages),
            modify_params_after_iter=adjust_lambda_a,
            iter_stop=2,
        )

    assert len(oct_result.info_vals) == 3
    assert isinstance(oct_result.info_vals[1], tuple)
    assert len(oct_result.info_vals[1]) == 2
    assert abs(oct_result.info_vals[1][0] - 0.001978333994757067) < 1e-8
    assert oct_result.info_vals[1][1] == 'λₐ: 0.5 → 0.25'
    assert 'Iteration 0: \tF = 0.000000' in stdout.getvalue()
    assert 'msg: λₐ: 1.0 → 0.5' in stdout.getvalue()
    assert 'Iteration 1: \tF = 0.001978' in stdout.getvalue()
    assert 'msg: λₐ: 0.5 → 0.25' in stdout.getvalue()


def test_invalid_print_table():
    J_T = lambda **kwargs: 1.0
    with pytest.raises(ValueError) as exc_info:
        print_table(J_T=J_T, col_formats=".2e")
    assert '8 elements' in str(exc_info.value)
    with pytest.raises(ValueError) as exc_info:
        print_table(J_T=J_T, col_formats="eeeeeeee")
    assert 'percent format string' in str(exc_info.value)
    with pytest.raises(ValueError) as exc_info:
        print_table(J_T=J_T, col_formats=(".2e", ".2e"))
    assert '8 elements' in str(exc_info.value)
    with pytest.raises(ValueError) as exc_info:
        # fmt: off
        col_formats = (
            '%d', '%.2', '%.2e', '%.2e', '%.2e', '%.2e', '%.2e', '%d',
        )
        # fmt: on
        print_table(J_T=J_T, col_formats=col_formats)
    assert 'Invalid col_formats' in str(exc_info.value)
    with pytest.raises(ValueError) as exc_info:
        print_table(J_T=J_T, col_headers="header")
    assert '8 elements' in str(exc_info.value)
    with pytest.raises(ValueError) as exc_info:
        print_table(J_T=J_T, col_headers=("header", "header"))
    assert '8 elements' in str(exc_info.value)

    # fmt: off
    col_headers = (
        "it", "J_T", "∫gₐ(ϵ{i})dt", "∑∫gₐ(t)dt", "J", "ΔJ_T", "ΔJ", "secs",
    )
    # fmt: on
    with pytest.raises(ValueError) as exc_info:
        print_table(
            J_T=J_T, show_g_a_int_per_pulse=True, col_headers=col_headers
        )
    assert "must support '.format(l=l)'" in str(exc_info.value)

    class no_format:
        """dummy broken g_a_hdr formatter"""

    # fmt: off
    col_headers = (
        "iteration", "J_T", no_format(), "∑∫gₐ(t)dt", "J", "ΔJ_T", "ΔJ",
        "secs",
    )
    # fmt: on
    with pytest.raises(ValueError) as exc_info:
        print_table(
            J_T=J_T, show_g_a_int_per_pulse=True, col_headers=col_headers
        )
    assert "must support '.format(l=l)'" in str(exc_info.value)


def test_print_table_custom_format(request, simple_state_to_state_system):
    """Test print_table with custom col_formats."""
    objectives, pulse_options, tlist = simple_state_to_state_system
    testdir = Path(os.path.splitext(request.module.__file__)[0])

    with io.StringIO() as log_fh:

        krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            store_all_pulses=True,
            info_hook=krotov.info_hooks.print_table(
                J_T=krotov.functionals.J_T_re,
                col_formats=(
                    '%8d',
                    '%12.4e',
                    '%12.4e',
                    '%12.4e',
                    '%12.4e',
                    '%12.4e',
                    '%12.4e',
                    '%05d',
                ),
                out=log_fh,
            ),
            iter_stop=1,
            skip_initial_forward_propagation=True,
        )

        log = log_fh.getvalue()

    expected_log = (testdir / 'custom_format_out.txt').read_text(
        encoding='utf8'
    )
    for (ln, ln_e) in zip(log.splitlines(), expected_log.splitlines()):
        # we have to remove the secs in the last column (not stable)
        assert ln[:-2] == ln_e[:-2]


def test_print_table_custom_header(request, simple_state_to_state_system):
    """Test print_table with custom col_formats."""
    objectives, pulse_options, tlist = simple_state_to_state_system
    testdir = Path(os.path.splitext(request.module.__file__)[0])

    with io.StringIO() as log_fh:

        krotov.optimize_pulses(
            objectives,
            pulse_options=pulse_options,
            tlist=tlist,
            propagator=krotov.propagators.expm,
            chi_constructor=krotov.functionals.chis_re,
            store_all_pulses=True,
            info_hook=krotov.info_hooks.print_table(
                J_T=krotov.functionals.J_T_re,
                col_headers=(
                    'iteration',
                    ' final time functional',
                    ' running cost (pulse {l})',
                    ' total running cost',
                    ' total functional',
                    ' change in final time functional',
                    ' change in total functional',
                    ' seconds for iteration',
                ),
                show_g_a_int_per_pulse=True,
                out=log_fh,
            ),
            iter_stop=1,
            skip_initial_forward_propagation=True,
        )

        log = log_fh.getvalue()

    expected_log = (testdir / 'custom_header_out.txt').read_text(
        encoding='utf8'
    )
    for (ln, ln_e) in zip(log.splitlines(), expected_log.splitlines()):
        # we have to remove the secs in the last column (not stable)
        assert ln[:-2] == ln_e[:-2]
