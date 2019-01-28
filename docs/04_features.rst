Features
========

* Simultaneously optimize over an arbitrary list of objectives
* Optimize over multiple control fields at the same time
* Arbitrary equations of motion, though a `propagator` callback function
* Arbitrary optimization functional, through `chi_constructor` callback function
* Allows injection of arbitrary code, through `modify_params_after_iter` function
* Customizable parallelization of the propagation of different objectives
* Customizable analysis and convergence check
* Support for dissipative dynamics (Liouville space)
* Convenience constructors for objectives describing gate optimization (in
  Hilbert- or Liouville space), and "ensemble optimization" for robust controls


**Not yet implemented:**

* non-linear controls
* state-dependent constraints
