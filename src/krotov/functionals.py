"""Functionals and `chi_constructor` routines"""

__all__ = ['chis_re']


def chis_re(states_T, objectives, tau_vals):
    r'''States $\ket{\chi}$ for the "real-part gate functional"'''
    c = float(1 / (2 * len(objectives)))
    return [c * obj.target_state for obj in objectives]
