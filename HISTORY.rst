=======
History
=======

1.2.1 (2021-01-13)
------------------

* Bugfix: Crash when initializing discretized numpy-array controls (`#79`_, thanks to `@loganbvh`_)
* Bugfix: Corrected definition of co-states in Dissipative Qubit Reset example (`#80`_, thanks to `Alberto Castro`_)
* Update: Switched Testing and Documentation deployment from Travis to Github Actions (`#82`_)

1.2.0 (2020-08-17)
------------------

* Added: ``via_midpoints`` argument to ``krotov.conversions.discretize`` function
* Changed: Controls and update shapes are now discretized in a way that ensures numerical stability (`#74`_, thanks to `@zachmanson`_)
* Changed: Replaced ``uniseg`` dependency with ``grapheme`` (`#76`_)

Note: due to the changes in the time discretization of the controls and update shapes, this version will generally not reproduce optimization results from previous versions to machine precision.


1.1.0 (2020-03-24)
------------------

* Added: Support for Python 3.8
* Added: Support for QuTiP 4.5.0
* Added: Support for parallelization with loky_ (`#72`_)
* Added: ``krotov.parallelization.set_parallelization`` function
* Added: ``krotov.parallelization.parallel_map`` function (improved implementation of QuTiP's ``parallel_map``)
* Added: Ability to use threadpoolctl_ to limit unwanted threading
* Added: `limit_thread_pool` option to ``krotov.optimize_pulses``
* Changed: ``krotov.propagators.expm`` now guarantees single-threaded execution


1.0.0 (2019-12-16)
------------------

* Update: Citation info now points to `SciPost paper <https://scipost.org/SciPostPhys.7.6.080>`_ (`#61`_)
* Added: parameters `col_formats` and `col_headers` to customize the output of ``krotov.info_hooks.print_table`` (`#65`_)
* Added: info-hooks now have access to the additional arguments `propagator`, `chi_constructor`, `mu`, `sigma`, `iter_start`, and `iter_stop` (`#66`_)
* Added: parameter `keep_original_objectives` to ``krotov.objectives.ensemble_objectives`` (`#67`_)
* Added: "Related Software" in the documentation
* Update: Documentation is now hosted on gh-pages_ and deployed by Doctr_ (`#68`_)


0.5.0 (2019-12-04)
------------------

* Update: Documentation now contains all information from https://arxiv.org/abs/1902.11284v5
* Added: Allow to pass `args` to time-dependent control functions (`#56`_, thanks to `@timohillmann`_)
* Changed: Renamed ``krotov.structural_conversions`` to ``krotov.conversions``
* Bugfix: Crash when ``krotov.optimize_pulses`` is called with ``iter_stop=0`` (`#58`_)
* Added: ``krotov.result.Result`` is now exposed at the top level of the API, as ``krotov.Result`` (`#59`_, thanks to `@nathanshammah`_)
* Added: str-representation of ``krotov.result.Result`` now includes the total running time (`#60`_, thanks to `@nathanshammah`_)


0.4.1 (2019-10-11)
------------------

* Update: Documentation now contains all information from https://arxiv.org/abs/1902.11284v4 (`#54`_)
* Added: a PDF of the documentation is now available at https://github.com/qucontrol/krotov/tree/master/docs/pdf (`#52`_, thanks to `@TejasAvinashShetty`_)


0.4.0 (2019-10-08)
------------------

* Added: Support for Python 3.7
* Changed: The ``'shape'`` key in ``pulse_options`` was renamed to ``'update_shape'``, to further avoid confusion between pulse shapes and update shapes.
* Changed: The ``.adjoint`` property of ``Objective`` is now a method
* Added: Ability to not use QuTiP ``Qobj`` objects, but arbitrary low-level objects instead.
* Improved: Printing an ``Objective`` now uses internal counters and a symbolic notation to identify objects shared between different objectives. (`#43`_)
* Improved: ``gate_objectives`` now takes into account if target states are (reshuffled) basis states and does not create unnecessary new copies.
* Bugfix: Two ``Objective`` instances that contain numpy arrays as controls can now be compared with ``==`` (`#44`_)
* Bugfix: Custom attributes (such as ``weight``) are now preserved when copying an ``Objective`` (`#44`_)
* Bugfix: Calling ``copy.deepcopy`` on an ``Objective`` now preserves control functions (`#44`_)
* Improved: The ``Objective.mesolve`` and ``Objective.propagate`` methods can now receive arguments ``H`` and ``c_ops`` to override the respective attributes of the objectives. This make is easier to analyze perform a robustness analysis, where the result of an optimization should be propagated under a perturbed Hamiltonian.
* Improved: The ``print_table`` and ``print_debug_information`` info-hooks now flush their output buffers after each iteration. As a result, when writing to a file, that file can be watched with ``tail -f``.
* Changed: Redefine ``tau_vals`` as their complex conjugate, fixing a bug in ``chis_ss`` and ``chis_sm`` (`#46`_)
* Bugfix: Correctly calculate ∂H/∂ϵ if ϵ occurs in H multiple times (`#47`_, thanks to `@uiofgh`_)
* Bugfix: Correctly calculate ∂H/∂ϵ=0 if the specific ϵ currently being updated does not occur in H (`#48`_)
* Added: Method ``objectives_with_controls`` for ``Result`` object.


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


.. _loky: https://loky.readthedocs.io/
.. _gh-pages: https://qucontrol.github.io/krotov
.. _Doctr: https://drdoctr.github.io
.. _threadpoolctl: https://github.com/joblib/threadpoolctl
.. _@uiofgh: https://github.com/uiofgh
.. _@TejasAvinashShetty: https://github.com/TejasAvinashShetty
.. _@timohillmann: https://github.com/timohillmann
.. _@nathanshammah: https://github.com/nathanshammah
.. _@zachmanson: https://github.com/zachmanson
.. _@loganbvh: https://github.com/loganbvh
.. _Alberto Castro: https://www.bifi.es/~acastro/
.. _#26: https://github.com/qucontrol/krotov/issues/26
.. _#29: https://github.com/qucontrol/krotov/issues/29
.. _#32: https://github.com/qucontrol/krotov/issues/32
.. _#35: https://github.com/qucontrol/krotov/issues/35
.. _#36: https://github.com/qucontrol/krotov/issues/36
.. _#39: https://github.com/qucontrol/krotov/issues/39
.. _#41: https://github.com/qucontrol/krotov/issues/41
.. _#43: https://github.com/qucontrol/krotov/issues/43
.. _#44: https://github.com/qucontrol/krotov/issues/44
.. _#46: https://github.com/qucontrol/krotov/issues/46
.. _#47: https://github.com/qucontrol/krotov/issues/47
.. _#48: https://github.com/qucontrol/krotov/issues/48
.. _#52: https://github.com/qucontrol/krotov/issues/42
.. _#54: https://github.com/qucontrol/krotov/issues/54
.. _#56: https://github.com/qucontrol/krotov/issues/56
.. _#58: https://github.com/qucontrol/krotov/issues/58
.. _#59: https://github.com/qucontrol/krotov/issues/59
.. _#60: https://github.com/qucontrol/krotov/issues/60
.. _#61: https://github.com/qucontrol/krotov/issues/61
.. _#65: https://github.com/qucontrol/krotov/issues/65
.. _#66: https://github.com/qucontrol/krotov/issues/66
.. _#67: https://github.com/qucontrol/krotov/issues/67
.. _#68: https://github.com/qucontrol/krotov/issues/68
.. _#72: https://github.com/qucontrol/krotov/issues/72
.. _#74: https://github.com/qucontrol/krotov/issues/74
.. _#76: https://github.com/qucontrol/krotov/issues/76
.. _#79: https://github.com/qucontrol/krotov/issues/79
.. _#80: https://github.com/qucontrol/krotov/issues/80
.. _#82: https://github.com/qucontrol/krotov/issues/82
