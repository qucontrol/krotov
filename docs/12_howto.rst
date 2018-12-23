How-Tos
=======

How to optimize towards a quantum gate
--------------------------------------

To optimize towards a quantum gate :math:`\Op{O}` in a *closed* quantum system,
set one :class:`.Objective` for state in the logical basis, with the basis
state :math:`\ket{\Psi_k}` as the :attr:`~.Objective.initial_state` and
:math:`\Op{O} \ket{\Psi_k}` as the :attr:`~.Objective.target`.

You may use :func:`krotov.gate_objectives() <krotov.objectives.gate_objectives>`
to construct the appropriate list of objectives.


How to optimize complex control fields
--------------------------------------

This implementation of Krotov's method requires real-valued control fields. You
must rewrite your Hamiltonian to contain the real part and the imaginary part
of the field as two independent controls. This is always possible. For example,
for a driven harmonic oscillator in the rotating wave approximation, the
interaction Hamiltonian is

.. math::

    \Op{H}_\text{int}
    = \epsilon^*(t) \Op{a} + \epsilon(t) \Op{a}^\dagger
    =  \epsilon_{\text{re}}(t) (\Op{a} + \Op{a}^\dagger) + \epsilon_{\text{im}}(t) (i \Op{a} - i \Op{a}^\dagger)

where :math:`\epsilon_{\text{re}}(t)= \Re[\epsilon(t)]` and
:math:`\epsilon_{\text{im}}(t) = \Im[\epsilon(t)]` are considered as two
independent (real-valued) controls.


How to exclude a control from the optimization
----------------------------------------------

In order to force the optimization to leave any particular control field
unchanged, set its update shape to
:func:`zero_shape <krotov.shapes.zero_shape>`::

    krotov.PulseOptions(lambda_a=1, shape=krotov.shapes.zero_shape)


How to define a new optimization functional
-------------------------------------------

In order to define a new optimization functional :math:`J_T`:

* Decide on what should go in :attr:`.Objective.target` to best describe the
  *physical* control target. If it the control target is fulfilled when the
  :attr:`.Objective.initial_state` evolves to a specific target state under the
  optimal control fields, that target state should to in
  :attr:`~.Objective.target`.

* Define a function `chi_constructor` that calculates the boundary
  condition for the backward-propagation in Krotov's method

  .. math::

        \ket{\chi_k(T)} \equiv - \left. \frac{\partial J_T}{\partial \bra{\phi_k(T)}} \right\vert_{\ket{\phi_k(T)}}

  or the equivalent experession in Liouville space. The function calculates the
  states :math:`\ket{\chi_k}` based  on the forward-propagated states
  $\ket{\phi_k(T)}$ and the list of objectives. For convenience, for when
  :attr:`~.Objective.target` contains a target state, `chi_constructor` will
  also receive `tau_vals` containing the overlaps
  :math:`\tau_k = \Braket{\phi_k(T)}{\phi_k^{\tgt}}`. See :func:`.chis_re` for
  an example.

* Optionally, define a function that can be used as an `info_hook`
  in :func:`.optimize_pulses` and that returns the value
  :math:`J_T`. This is not required to run an optimization, as the
  functional is entirely implicit in `chi_constructor`. However, calculating
  the value of the functional is useful for convergence analysis
  (`check_convergence` in :func:`.optimize_pulses`)


How to optimize towards a two-qubit gate up to single-qubit corrections
-----------------------------------------------------------------------

Use :func:`krotov.objectives.gate_objectives` with ``local_invariants=True`` in
order to construct a list of objectives suitable for an optimization using a
"local-invariant functional" :cite:`MullerPRA11`.

The |weylchamber package|_ contains the suitable `chi_constructor` routines to
pass to :func:`.optimize_pulses`.



How to penalize population in a forbidden subspace
--------------------------------------------------

In principal, :func:`.optimize_pulses` has a `state_dependent_constraint`.
However, this has some caveats. Most notably, it results in an inhomogeneous
equation of motion, which is currently not implemented.

The recommended "workaround" is to place artificially high dissipation on the
the levels in the forbidden subspace. A non-Hermitian Hamiltonian is usually a
good way to realize this.


How to optimize towards an arbitrary perfect entangler
------------------------------------------------------

Use :func:`krotov.objectives.gate_objectives` with ``gate=PE`` in
order to construct a list of objectives suitable for an optimization using a
"perfect entanglers" :cite:`WattsPRA2015,GoerzPRA2015`.

The |weylchamber package|_ contains the suitable `chi_constructor` routines to
pass to :func:`.optimize_pulses`.

.. |weylchamber package| replace:: ``weylchamber`` package
.. _weylchamber package: https://github.com/qucontrol/weylchamber
