"""Functionals and `chi_constructor` routines"""

__all__ = ['chis_re']


def chis_re(states_T, objectives, tau_vals):
    r'''States $\ket{\chi}$ for the "real-part gate functional"'''
    c = float(1 / (2 * len(objectives)))
    # `obj.target` is assumed to be the "target state" (gate applied to
    # `initial_state`)
    return [c * obj.target for obj in objectives]
