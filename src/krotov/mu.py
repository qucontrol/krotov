r"""Routines for `mu` in :func:`krotov.optimize.optimize_pulses`.

The first-order Krotov update equation is usually written as

.. math::

       \Delta \epsilon(t) \propto \Im
       \Bigg\langle
         \chi_k^{(i)}(t)
       \Bigg\vert
         \Bigg(
            \left.\frac{\partial \Op{H}}{\partial \epsilon}\right\vert_{
                {\scriptsize \begin{matrix}
                    \phi^{(i+1)}(t)\\\epsilon^{(i+1)}(t)
                \end{matrix}}}
         \Bigg)
       \Bigg\vert
         \phi_k^{(i+1)}(t)
       \Bigg\rangle\,,

where $\ket{\chi_k}$ are states backward-propagated from a boundary condition
determined by the functional, $\ket{\phi_k}$ are forward-propagated from the
initial states, and $\frac{\partial \Op{H}}{\partial \epsilon}$ is the
derivative of the Hamiltonian with respect to the field. However, this is true
only for Hilbert-space states evolving under a standard Schrödinger equation.

More generally (e.g. when the states $\chi_k$ and $\phi_k$ are density matrices
and the equation of motion is the master equation in Lindblad form), the
correct formulation is

.. math::

    \frac{\partial \Op{H}}{\partial \epsilon}
    \rightarrow
    \mu = \frac{\partial H}{\partial \epsilon}\,,

where $H$ is now the abstract operator appearing in the equation of motion of
the abstract state

.. math::

    \dot\phi_k(t) = -i H \phi_k(t)


For density matrices (using the most common convention for the Liouvillian
$\Liouville$, and the one that QuTiP makes in
:func:`qutip.superoperator.liouvillian`), we have

.. math::

    \frac{\partial}{\partial t}\Op{\rho}_k(t) = \Liouville \Op{\rho}_k(t)

and thus $H = i \Liouville$.

To allow for arbitrary equations of motion, a routine `mu` may be passed to
:func:`.optimize_pulses` that returns the abstract operator $\mu$ as a
:class:`~qutip.Qobj`, or alternatively as a callable that takes $\phi_k$ as its
argument and evaluates $\mu \phi_k$. The default `mu` is
:func:`derivative_wrt_pulse`, which associates the abstract `H` directly with
the attribute `H` of each `objective`. Alternative implementations of `mu` must
have the same signature.
"""

__all__ = ['derivative_wrt_pulse']


def derivative_wrt_pulse(objective, pulses, mapping, i_pulse):
    r"""Calculate ∂H/∂ϵ directly from the `H` attribute of the `objective`"""
    ham_mapping = mapping[0][i_pulse]
    if len(ham_mapping) == 0:
        return 0
    else:
        mu = objective.H[ham_mapping[0]][0]
        for i in ham_mapping[1:]:
            mu += objective.H[ham_mapping[i]][0]
    for i_c_op in range(len(objective.c_ops)):
        if len(mapping[i_c_op+1][i_pulse]) != 0:
            raise NotImplementedError(
                "Time-dependent collapse operators not implemented")
    return mu
