=======
History
=======

(next version)
--------------

* Improved: Printing an ``Objective`` now uses internal counters and a symbolic notation to identify objects shared between different objectives. (`#43`_)
* Improved: ``gate_objectives`` now takes into account if target states are (reshuffled) basis states and does not create unnecessary new copies.
* Bugfix: Two ``Objective`` instances that contain numpy arrays as controls can now be compared with ``==`` (`#44`_)
* Bugfix: Custom attributes (such as ``weight``) are now preserverd when copying an ``Objective`` (`#44`_)
* Bugfix: Calling ``copy.deepcopy`` on an ``Objective`` now preserves control functions (`#44`_)
* Improved: The ``Objective.mesolve`` and ``Objective.propagate`` methods can now receive arguments ``H`` and ``c_ops`` to override the respective attributes of the objectives. This make is easier to analyze perform a robustness analysis, where the result of an optimization should be propagated under a perturbed Hamiltonian.
* Improved: The ``print_table`` and ``print_debug_information`` info-hooks now flush their output buffers after each iteration. As a result, when writing to a file, that file can be watched with ``tail -f``.

0.3.0 (2019-03-01)
------------------

* Added: Preprint citation information (``krotov.__arxiv__``, ``krotov.__citation__``, ``krotov.__bibtex__``)
* Added: Ability to continue from a previous optimization (`#26`_)
* Added: Parameter ``out`` to ``print_table`` info-hook
* Added: Parameter ``finalize`` to ``Result.load``
* Added: Ability to dump optimization result every so many iterations (``dump_result`` check-convergence routine)
* Added: `re-entrant` option for ``DensityMatrixODEPropagator``
* Bugfix: Discretize controls to float values (`#41`_)
* Bugfix: Fix overlap for non-Hermitian operators (`#39`_)
* Bugfix: Interface for passing ``tau_vals`` to ``chi_constructor`` (`#36`_)
* Added: function ``above_value`` for convergence check (`#35`_)


0.2.0 (2019-02-14)
------------------

* Added: Implementation of all the standard functionals
* Added: The ``info_hook`` receives additional information, including ∫gₐ(t)dt (`#32`_)
* Added: Initialization of objectives for gate optimization in Liouville space
* Added: A new propagator ``DensityMatrixODEPropagator`` for faster density matrix propagation
* Added: Support for "stateful" propagators by subclassing from ``krotov.propagators.Propagator``
* Changed: more flexibility for parallelization (`#29`_)
* Added: Support for the second-order pulse update
* Changed: The options for the controls (λₐ, update-shape) are now passed through a simplified ``dict`` interface, instead of a custom ``PulseOptions`` class.


0.1.0 (2018-12-24)
------------------

* Initial release with complete implementation of first-order Krotov's method
* Support for state-to-state and gate optimization, for both closed and open systems


.. _#26: https://github.com/qucontrol/krotov/issues/26
.. _#29: https://github.com/qucontrol/krotov/issues/29
.. _#32: https://github.com/qucontrol/krotov/issues/32
.. _#35: https://github.com/qucontrol/krotov/issues/35
.. _#36: https://github.com/qucontrol/krotov/issues/36
.. _#39: https://github.com/qucontrol/krotov/issues/39
.. _#41: https://github.com/qucontrol/krotov/issues/41
.. _#43: https://github.com/qucontrol/krotov/issues/43
.. _#44: https://github.com/qucontrol/krotov/issues/44
