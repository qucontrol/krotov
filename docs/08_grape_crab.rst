Krotov vs GRAPE and CRAB
========================

Gradient-based optimization
---------------------------

At its core, Krotov's method is a gradient-based optimization method, and
most directly compares to
GRAdient-based-Pulse-Engineering (GRAPE), another gradient-based method
widely used in quantum control :cite:`KhanejaJMR05`. Historically, Krotov's
method has been widely used in the control of atomic and molecular dynamics,
whereas GRAPE has its origin in NMR. Generally, both methods can be used
interchangeably, see :ref:`choosing-an-optimization-method` below for some of
the subtleties.

GRAPE is included in QuTiP, see the `section on Quantum Optimal Control in the QuTiP docs`_.
It is used via the :func:`qutip.control.pulseoptim.optimize_pulse` function.

This function has some limitations compared to the :func:`.optimize_pulses`
provided by the Krotov package. Most importantly:

* QuTiP's GRAPE implementation only allows for a single control field (hence *"pulse"*, instead of *"pulses"* in :func:`.optimize_pulses`)
* The optimization for multiple simultaneous objectives in GRAPE is limited to optimizing a quantum gate. Other uses of parallel objectives, such as optimizing for robustness, is not available.
* The GRAPE implementation does not allow for an arbitrary guess-pulse. Guess pulses can only be chosen from a specific set of options (including "random")
* There is only a fixed number of choices for choosing a method for time-propagation. Supplying a problem-specific propagator is not possible.
* The optimized pulses are defined on the intervals of the time grid, which does not fit in with other parts of QuTiP, specifically :func:`qutip.mesolve.mesolve`.

All of these are restrictions of the implementation, not of the underlying method.

.. note::
    There are plans to provide a wrapper-function that provides the interface
    of :func:`qutip.control.pulseoptim.optimize_pulse` for use with Krotov's
    method.

.. _section on Quantum Optimal Control in the QuTiP docs: http://qutip.org/docs/latest/guide/guide-control.html

Gradient-free optimization
--------------------------

In situations where the controls can be reduced to a relatively small number of
controllable parameters, gradient-free optimization becomes feasible. Gradient-free
optimization requires at most half of the numerical effort than either Krotov's
method or GRAPE, as they only require a single forward-propagation of the
objectives in each iteration. Such a gradient-free optimization can be easily realized
by following two simple steps:

* Write a function that takes an array of optimization parameters as input
  and returns a figure of merit on output. This function would, e.g., construct a
  numerical control pulse from the control parameters, simulate the dynamics
  using :func:`qutip.mesolve.mesolve`, and compare to a target state by evaluating
  the state overlap as a figure of merit.
* Pass the function to :func:`scipy.optimize.minimize` for gradient-free optimization.

The implementation in :func:`scipy.optimize.minimize` allows to choose between
different optimization methods, with "Nelder-Mead simplex" being the the default.

A special case of gradient-free optimization is the Chopped RAndom Basis (CRAB)
method :cite:`DoriaPRL11,CanevaPRA2011`.
The essence of CRAB is in the specific choice of the parametrization in terms of a
low-dimensional random basis, as the name implies. The optimization itself is performed
by Nelder-Mead simplex using this parametrization.

An implementation of CRAB is included in QuTiP, see `QuTiP's documentation of
CRAB`_, and uses the same :func:`qutip.control.pulseoptim.optimize_pulse`
interface as the GRAPE method discussed above with the same limitations.

.. _QuTiP's documentation of CRAB: http://qutip.org/docs/latest/guide/guide-control.html#the-crab-algorithm


.. _choosing-an-optimization-method:

Choosing an optimization method
-------------------------------

Whether to use a gradient-free optimization method, GRAPE, or
Krotov's method depends on the size of the problem (both the Hilbert
space dimension and the number of control parameters), the requirements
on the control pulse, and the optimization functional. Gradient-free
methods should be used if propagation is extremely cheap (small Hilbert
space dimension), the number of independent control parameters is
relatively small, or the functional is of a form that does not allow to
calculate gradients easily.

GRAPE should be used if the control parameters are discrete,
such as on a coarse-grained time grid, and the derivative of :math:`J`
with respect to each control parameter is easily computable. Moreover, evaluation of
the gradient should be ideally numerically inexpensive.

Krotov's method should be used if the control is near-continuous, and if
the derivative of :math:`J_T` with respect to the states, Eq. :eq:`chi_boundary`, can be
easily calculated. When these conditions are met, Krotov's method gives excellent convergence,
although it is often observed to slow down when getting close to the
minimum of :math:`J`. It can be beneficial to switch from Krotov's method to
GRAPE with LBFGS-B in the final stage of the optimization which does not show
such a slow-down due to inclusion of second-derivative information by the
Hessian on top of the gradient.

.. .. bibliography:: refs.bib
   :cited:
   :style: unsrt
