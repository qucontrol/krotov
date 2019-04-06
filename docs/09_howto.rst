How-Tos
=======

How to optimize towards a quantum gate
--------------------------------------

To optimize towards a quantum gate :math:`\Op{O}` in a *closed* quantum system,
set one :class:`.Objective` for each state in the logical basis, with the basis
state :math:`\ket{\phi_k}` as the :attr:`~.Objective.initial_state` and
:math:`\Op{O} \ket{\phi_k}` as the :attr:`~.Objective.target`.

You may use :func:`krotov.gate_objectives() <krotov.objectives.gate_objectives>`
to construct the appropriate list of objectives. See the
:ref:`/notebooks/05_example_transmon_xgate.ipynb` for an example. For more
advanced gate optimizations, also see :ref:`HowtoLIOptimization`,
:ref:`HowtoPEOptimization`, :ref:`HowtoDissipativeOptimization`, and
:ref:`HowtoRobustOptimization`.


How to optimize complex control fields
--------------------------------------

This implementation of Krotov's method requires real-valued control fields. You
must rewrite your Hamiltonian to contain the real part and the imaginary part
of the field as two independent controls. This is always possible. For example,
for a driven harmonic oscillator in the rotating wave approximation, the
interaction Hamiltonian is given by

.. math::

    \Op{H}_\text{int}
    = \epsilon^*(t) \Op{a} + \epsilon(t) \Op{a}^\dagger
    =  \epsilon_{\text{re}}(t) (\Op{a} + \Op{a}^\dagger) + \epsilon_{\text{im}}(t) (i \Op{a}^\dagger - i \Op{a})\,,

where :math:`\epsilon_{\text{re}}(t)= \Re[\epsilon(t)]` and
:math:`\epsilon_{\text{im}}(t) = \Im[\epsilon(t)]` are considered as two
independent (real-valued) controls.

See the :ref:`/notebooks/02_example_lambda_system_rwa_complex_pulse.ipynb` for an example.


How to stop the optimization when the error crosses some threshold
------------------------------------------------------------------

By default, an optimization stops after a predefined number of iterations
(`iter_stop` parameter in :func:`.optimize_pulses`). However, through the
interplay of the `info_hook` and the `check_convergence` routine  passed to
:func:`.optimize_pulses`, the optimization can be stopped based on the
optimization success or the rate of convergence: The `info_hook` routine should
return the value of the optimization functional or error, which is accessible to
`check_convergence` via the :attr:`.Result.info_vals` attribute, see
:mod:`krotov.convergence` for details.

Generally, you should use the :func:`krotov.info_hooks.print_table` function as
an `info_hook`, which receives a function to evaluate the optimization
functional $J_T$ as a parameter. Then, use
:func:`krotov.convergence.value_below` as a `check_convergence` routine to stop
the optimization when $J_T$ falls below some given threshold.

See the :ref:`/notebooks/02_example_lambda_system_rwa_complex_pulse.ipynb` for
an example.


How to exclude a control from the optimization
----------------------------------------------

In order to force the optimization to leave any particular control field
unchanged, set its update shape to :func:`krotov.shapes.zero_shape`
in the `pulse_options` that you pass to :func:`.optimize_pulses`.


How to define a new optimization functional
-------------------------------------------

In order to define a new optimization functional :math:`J_T`:

* Decide on what should go in :attr:`.Objective.target` to best describe the
  *physical* control target. If the control target is reached when the
  :attr:`.Objective.initial_state` evolves to a specific target state under the
  optimal control fields, that target state should be included in
  :attr:`~.Objective.target`.

* Define a function `chi_constructor` that calculates the boundary
  condition for the backward-propagation in Krotov's method,

  .. math::

        \ket{\chi_k(T)} \equiv - \left. \frac{\partial J_T}{\partial \bra{\phi_k(T)}} \right\vert_{\ket{\phi_k(T)}}\,,

  or the equivalent experession in Liouville space. This function should
  calculate the states :math:`\ket{\chi_k}` based  on the forward-propagated
  states :math:`\ket{\phi_k(T)}` and the list of objectives. For convenience,
  when :attr:`~.Objective.target` contains a target state, `chi_constructor`
  will also receive `tau_vals` containing the overlaps :math:`\tau_k =
  \Braket{\phi_k(T)}{\phi_k^{\tgt}}`. See :func:`.chis_re` for an example.

* Optionally, define a function that can be used as an `info_hook`
  in :func:`.optimize_pulses` which returns the value
  :math:`J_T`. This is not required to run an optimization since the
  functional is entirely implicit in `chi_constructor`. However, calculating
  the value of the functional is useful for convergence analysis
  (`check_convergence` in :func:`.optimize_pulses`)

See :mod:`krotov.functionals` for some standard functionals. An example for a
more advanced functional is the :ref:`/notebooks/07_example_PE.ipynb`.


How to penalize population in a forbidden subspace
--------------------------------------------------

In principle, :func:`.optimize_pulses` has a `state_dependent_constraint`.
However, this has some caveats. Most notably, it results in an inhomogeneous
equation of motion, which is currently not implemented.

The recommended "workaround" is to place artificially high dissipation on the
levels in the forbidden subspace. A non-Hermitian Hamiltonian is usually a
good way to realize this. See the
:ref:`/notebooks/03_example_lambda_system_rwa_non_hermitian.ipynb`
for an example.


.. _HowtoLIOptimization:

How to optimize towards a two-qubit gate up to single-qubit corrections
-----------------------------------------------------------------------

Use :func:`krotov.objectives.gate_objectives` with ``local_invariants=True`` in
order to construct a list of objectives suitable for an optimization using a
"local-invariant functional" :cite:`MullerPRA11`. This optimizes towards a
point in the `Weyl chamber`_.

The |weylchamber package|_ contains the suitable `chi_constructor` routines to
pass to :func:`.optimize_pulses`.


.. _HowtoPEOptimization:

How to optimize towards an arbitrary perfect entangler
------------------------------------------------------

Closely related to an optimization towards a point in the Weyl chamber is the
optimization towards an arbitrary perfectly entangling two-qubit gate.
Geometrically, this means optimizing towards the polyhedron of perfect
entanglers in the `Weyl chamber`_.

Use :func:`krotov.objectives.gate_objectives` with ``gate='PE'`` in
order to construct a list of objectives suitable for an optimization using a
"perfect entanglers" functional :cite:`WattsPRA2015,GoerzPRA2015`.
This is illustrated in the :ref:`/notebooks/07_example_PE.ipynb`.

Again, the `chi_constructor` is available in the |weylchamber package|_.

.. |weylchamber package| replace:: ``weylchamber`` package
.. _weylchamber package: https://github.com/qucontrol/weylchamber
.. _Weyl chamber: https://weylchamber.readthedocs.io/en/latest/tutorial.html


.. _HowtoDissipativeOptimization:

How to optimize in a dissipative system
---------------------------------------

To optimize a dissipative system, it is sufficient to set an :class:`.Objective`
with a density matrix for the :attr:`~.Objective.initial_state` and
:attr:`~.Objective.target`, and a Liouvillian in :attr:`.Objective.H`.
See the :ref:`/notebooks/04_example_dissipative_qubit_reset.ipynb` for an
example.

Instead of a Liouvillian, it is also possible to set :attr:`.Objective.H` to
the system Hamiltonian, and :attr:`.Objective.c_ops` to the appropriate
Lindblad operators. However, it is generally much more efficient to use
:func:`krotov.objectives.liouvillian` to convert a time-dependent Hamiltonian
and a list of Lindblad operators into a time-dependent Liouvillian. In either
case, the `propagate` routine passed to :func:`~krotov.optimize.optimize_pulses`
must be aware of and compatible with the convention for the objectives.

Specifically for gate optimization, the routine
:func:`~krotov.objectives.gate_objectives`
can be used to automatically set appropriate objectives for an optimization in
Liouville space. The parameter `liouville_states_set` indicates that the system
dynamics are in Liouville space and sets an appropriate choice of matrices that
track the optimization according to Ref. :cite:`GoerzNJP2014`.
See the :ref:`/notebooks/06_example_3states.ipynb` for an example.

For weak dissipation, it may also be possible to avoid the use of density
matrices altogether, and to instead use a non-Hermitian Hamiltonian. For example, you may
use the effective Hamiltonian from the MCWF method :cite:`PlenioRMP1998`,

.. math::

   \Op{H}_{\text{eff}} = \Op{H} - \frac{i}{2} \sum_k \Op{L}_k^\dagger \Op{L}_k\,,

for the Hermitian Hamiltonian :math:`\Op{H}` and the Lindblad operators
:math:`\Op{L}_k`.  Propagating :math:`\Op{H}_{\text{eff}}` (without quantum
jumps) will lead to a decay in the norm of the state corresponding to how much
dissipation the state is subjected to. Numerically, this will usually increase
the value of the optimization functional (that is, the error). Thus the
optimization can be pushed towards avoiding decoherence, without explicitly
performing the optimization in Liouville space. See the
:ref:`/notebooks/03_example_lambda_system_rwa_non_hermitian.ipynb` for an
example.


.. _HowtoRobustOptimization:

How to optimize for robust pulses
---------------------------------

Control pulses can be made robust with respect to variations in the system by
doing an ensemble optimization, as proposed in Ref. :cite:`GoerzPRA2014`. The
idea if to sample a representative selection of possible system Hamiltonians,
and to optimize over an *average* of the entire ensemble.

An appropriate set of objectives can be generated with the
:func:`~krotov.objectives.ensemble_objectives` function.


How to parallelize the optimization
-----------------------------------

Krotov's method is inherently parallel across different objectives. See
:mod:`krotov.parallelization`, and the
:ref:`/notebooks/05_example_transmon_xgate.ipynb` for an example.

.. _HowtoStoreResult:

How to prevent losing an optimization result
--------------------------------------------

Optimizations usually take several hundred to several thousand iterations to
fully converge. Thuse, the :func:`.optimize_pulses` routine  may require
significant runtime (often multiple days for large problems). Once an
optimization has completed, you are strongly encouraged to store the result to
disk, using :meth:`.Result.dump`.  You may also consider using
:func:`.dump_result` during the `check_convergence` step to dump the current
state of the optimization to disk at regular intervals. This protects you from
losing work if the optimization is interrupted in any way, like an unexpected
crash.

In order to continue after such a crash, you can restore a :class:`.Result`
object containing the recent state of the optimization using
:meth:`.Result.load` (with the original `objectives` and ``finalize=True`` if
the dump file originates from :func:`.dump_result`). You may then call
:func:`.optimize_pulses` and pass the loaded :class:`.Result` object as
`continue_from`.  The new optimization will start from the most recent
optimized controls as a guess, and continue to count iterations from the
previous result. See :ref:`HowtoContinueOptimization` for further details.


.. _HowtoContinueOptimization:

How to continue from a previous optimization
--------------------------------------------

See :ref:`HowtoStoreResult` for how to continue from an optimization that ended
(crashed) prematurely.  Even when an optimization has completed normally, you
may still want to continue with further iterations -- either because you find
that the original `iter_stop` was insufficient to reach full convergence, or
because you would like to modify some parameters, like the λₐ values for
each control. In this case, you can again call :func:`.optimize_pulses` and
pass the :class:`.Result` object from the previous optimization as
`continue_from`. Note that while you are free to change the `pulse_options`
between the two optimization, the `objectives` must remain the same. The
functional (`chi_constructor`) and the `info_hook` should also remain the same
(otherwise, you may and up with inconsistencies in your :class:`.Result`). The
:class:`.Result` object returned by the second optimization will include all
the data from the first optimization.


How to maximize numerical efficiency
------------------------------------

For systems of non-trivial size, the main numerical effort should be in the
simulation of the system dynamics. Every iteration of Krotov's method requires
a full backward propagation and a full forward propagation of the states associated with each
objective, see :mod:`krotov.propagators`. Therefore, the best numerical
efficiency can be achieved by optimizing the performance of the `propagator`
that is passed to :func:`~krotov.optimize.optimize_pulses`.

One possibility is to implement problem-specific propagators, such as
:class:`krotov.propagators.DensityMatrixODEPropagator`. Going further, you
might consider implementing the propagator with the help of lower-level instructions, e.g.,
by using Cython_.

.. _Cython: https://cython.org


How to deal with the optimization running out of memory
-------------------------------------------------------

Krotov's method requires the storage of at least one set of propagated state
over the entire time grid, for each objective. For the second-order update
equation, up to three sets of stored states per objective may be required. In
particular for larger systems and dynamics in Liouville space, the memory
required for storing these states may be prohibitively expensive.

The :func:`~krotov.optimize.optimize_pulses` accepts a `storage` parameter
to which a constructor for an array-like container can be passed wherein the
propagated states will be stored. It is possible to pass custom out-of-memory
storage objects, such as Dask_ arrays. This may carry a significant penalty in
runtime, however, as states will have to be read from disk, or across the
network.

.. _Dask: http://docs.dask.org/en/latest/
