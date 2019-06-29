"""Test for modify_params_after_iter/info_hook"""
import contextlib
import io

import numpy as np

import krotov
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
