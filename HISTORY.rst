=======
History
=======


(next release)
--------------

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


.. _#29: https://github.com/qucontrol/krotov/issues/29
.. _#32: https://github.com/qucontrol/krotov/issues/32
.. _#35: https://github.com/qucontrol/krotov/issues/35
.. _#36: https://github.com/qucontrol/krotov/issues/36
.. _#39: https://github.com/qucontrol/krotov/issues/39
.. _#41: https://github.com/qucontrol/krotov/issues/41
