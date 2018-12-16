How-Tos
=======

How to optimize towards a quantum gate
--------------------------------------

To optimize towards a quantum gate :math:`\Op{O}` in a *closed* quantum system,
set one :class:`.Objective` for state in the logical basis, with the basis
state :math:`\ket{\Psi_k}` as the `initial_state` and :math:`\Op{O}
\ket{\Psi_k}` as the `target_state`.

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


How to define an optimization functional
----------------------------------------

TODO


How to optimize towards a two-qubit gate up to single-qubit corrections
-----------------------------------------------------------------------

TODO


How to penalize population in a forbidden subspace
--------------------------------------------------

TODO


How to optimize towards an arbitrary perfect entangler
------------------------------------------------------

TODO


How to optimize in a dissipative system
---------------------------------------

TODO


How to optimize for robust pulses
---------------------------------

TODO


How to parallelize the optimization
-----------------------------------

TODO


How to maximize numerical efficiency
------------------------------------

TODO


How to deal with the optimization running out of memory
-------------------------------------------------------

TODO
