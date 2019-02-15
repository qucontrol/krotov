=======
History
=======


(next release)
--------------

* Added: function ``above_value`` for convergence check (#35)


0.2.0 (2019-02-14)
------------------

* Added: Implementation of all the standard functionals
* Added: The ``info_hook`` receives additional information, including ∫gₐ(t)dt (#32)
* Added: Initialization of objectives for gate optimization in Liouville space
* Added: A new propagator ``DensityMatrixODEPropagator`` for faster density matrix propagation
* Added: Support for "stateful" propagators by sublassing from ``krotov.propagators.Propagator``
* Changed: more flexibility for parallelization (#29)
* Added: Support for the second-order pulse update
* Changed: The options for the controls (λₐ, update-shape) are now passed through a simplified ``dict`` interface, instead of a custom ``PulseOptions`` class.


0.1.0 (2018-12-24)
------------------

* Initial release with complete implementation of first-order Krotov's method
* Support for state-to-state and gate optimization, for both closed and open systems


0.0.1 (2018-11-06)
------------------

* Non-functional placeholder release
