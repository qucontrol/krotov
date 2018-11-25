"""Routines that can be passed as `propagator` to :func:`.optimize_pulses`"""

import qutip

__all__ = ['expm']


def _apply(A, b):
    """Apply abstract operator `A` to state `b`."""
    if not(isinstance(A, qutip.Qobj)) and callable(A):
        return A(b)
    elif A.type == 'super' and b.type == 'oper':
        # Workaround for https://github.com/qutip/qutip/issues/939
        return qutip.vector_to_operator(A * qutip.operator_to_vector(b))
    else:
        return A * b


def expm(H, state, dt, c_ops, backwards=False):
    """Propagate using matrix exponentiation"""
    if len(c_ops) > 0:
        raise NotImplementedError("Liouville exponentiation not implemented")
    assert isinstance(H, list) and len(H) > 0
    if backwards:
        dt = -dt
    eqm_factor = -1j  # factor in front of H on rhs of the equation of motion
    if isinstance(H[0], list):
        if H[0][1].type == 'super':
            eqm_factor = 1
        A = (eqm_factor * H[0][1]) * H[0][0]
    else:
        if H[0].type == 'super':
            eqm_factor = 1
        A = eqm_factor * H[0]
    for part in H[1:]:
        if isinstance(part, list):
            A += (eqm_factor * part[1]) * part[0]
        else:
            A += (eqm_factor * part)
    ok_types = (
        (state.type == 'oper' and A.type == 'super') or
        (state.type in ['ket', 'bra'] and A.type == 'oper'))
    if ok_types:
        return _apply((A * dt).expm(), state)
    else:
        raise NotImplementedError(
            "Cannot handle argument types A:%s, state:%s"
            % (A.type, state.type))
