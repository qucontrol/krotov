"""Routines that can be passed as `propagator` to :func:`.optimize_pulses`"""

__all__ = ['propagator_expm']


def propagator_expm(H, state, dt, c_ops):
    """Propagate using matrix exponentiation"""
    if len(c_ops) > 0:
        raise NotImplementedError("Liouville exponentiation not implemented")
    if not state.isket:
        raise NotImplementedError("Propagation of non-kets not implemented")
    assert isinstance(H, list) and len(H) > 0
    if isinstance(H[0], list):
        A = (-1j * H[0][1]) * H[0][0]
    else:
        A = -1j * H[0]
    for part in H[1:]:
        if isinstance(part, list):
            A += (-1j * part[1]) * part[0]
        else:
            A += (-1j * part)
    return (A * dt).expm() * state
