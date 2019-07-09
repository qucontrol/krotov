r"""Routines for `mu` in :func:`krotov.optimize.optimize_pulses`

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
only for Hilbert-space states evolving under a Schrödinger equation.

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


For density matrices, we have

.. math::

    \frac{\partial}{\partial t}\Op{\rho}_k(t) = \Liouville \Op{\rho}_k(t)

and thus $H = i \Liouville$.

To allow for arbitrary equations of motion, a routine `mu` may be passed to
:func:`.optimize_pulses` that returns the abstract operator $\mu$ as a
:class:`~qutip.Qobj`, or alternatively as a callable that takes $\phi_k$ as its
argument and evaluates $\mu \phi_k$. The default `mu` is
:func:`derivative_wrt_pulse`, which covers the most common equation of motions:

* standard Schrödinger equation
* master equation, where either the `H` attribute of the objective contains a
  Hamiltonian and there are Lindblad operators in `c_ops`, or the `H` attribute
  contains a super-operator $\Liouville$ directly (the case discussed above).

Alternative implementations of `mu` must have the same signature as
:func:`derivative_wrt_pulse`, but should only be required in rare
circumstances, such as when the derivative still depends on the control values
or on the states. (Or, if you can provide a more efficient problem-specific
implementation).
"""

__all__ = ['derivative_wrt_pulse']


def derivative_wrt_pulse(
    objectives, i_objective, pulses, pulses_mapping, i_pulse, time_index
):
    r"""Calculate ∂H/∂ϵ for the standard equations of motion.

    Args:
        objectives (list): List of :class:`.Objective` instances
        i_objective (int): The index of the objective in `objectives` whose
            equation of motion the derivative should be calculated.
        pulses (list): The list of pulses occuring in `objectives`
        pulses_mapping (list): The mapping of elements of `pulses` to the
            components of `objectives`, as returned by
            :func:`.extract_controls_mapping`
        i_pulse (int): The index of the pulse in `pulses` for which to
            calculate the derivative
        time_index (int): The index of the value in ``pulses[i_pulse]`` that
            should be plugged in to ∂H/∂ϵ. Not used, as this routine only
            considers equations of motion that are linear in the controls.

    Returns:
        callable: The quantum operator or super-operator that
        represents ∂H/∂ϵ. In general, the return type can be any callable `mu`
        so that ``mu(state)`` calculates the result of applying ∂H/∂ϵ to
        `state`. In most cases, a :class:`~qutip.Qobj` will be returned, which
        is just the most convenient example of an appropriate callable.

    This function covers the following cases:

    * the :attr:`~.Objective.H` attribute of the objective contains a
      Hamiltonian, there are no :attr:`~.Objective.c_ops` (Schrödinger
      equation: the abstract H in ∂H/∂ϵ is the Hamiltonian directly)

    * the :attr:`~.Objective.H` attribute of the objective contains a
      Hamiltonian $\Op{H}$, and there are Lindblad operators $\Op{L}_i$ in
      :attr:`~.Objective.c_ops` (master equation in Lindblad form). The
      abstract H is $i \Liouville$ for the Liouvillian defined as

      .. math::

        \Liouville[\Op{\rho}] =
        -i[\Op{H},\Op{\rho}]+\sum_{i} \left(
            \Op{L}_i \Op{\rho} \Op{L}_i^\dagger -
            \frac{1}{2} \left\{
                \Op{L}_i^\dagger \Op{L}_i, \Op{\rho}\right\} \right)

    * the :attr:`~.Objective.H` attribute of the objective contains a
      super-operator $\Liouville$, there are no :attr:`~.Objective.c_ops`
      (general master equation). The abstract H is again $i \Liouville$.
    """
    objective = objectives[i_objective]
    ham_mapping = pulses_mapping[i_objective][0][i_pulse]
    eqm_factor = -1j  # the factor in front of objective.H in the eqm
    if len(ham_mapping) == 0:
        return lambda state: 0 * state
    else:
        mu = objective.H[ham_mapping[0]][0]
        if mu.type == 'super':
            eqm_factor = 1
            mu *= 1j
        for i in ham_mapping[1:]:
            mu += (1j * eqm_factor) * objective.H[i][0]
    for i_c_op in range(len(objective.c_ops)):
        if len(pulses_mapping[i_objective][i_c_op + 1][i_pulse]) != 0:
            raise NotImplementedError(
                "Time-dependent collapse operators not implemented"
            )
    return mu
