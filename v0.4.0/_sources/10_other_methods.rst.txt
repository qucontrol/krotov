Other Optimization Methods
==========================

At its core, Krotov’s method is a gradient-based optimization method,
and most directly compares to GRadient Ascent Pulse Engineering
(GRAPE) :cite:`KhanejaJMR05`, another gradient-based method
widely used in quantum control :cite:`GlaserEPJD2015`. We
therefore compare to GRAPE first and highlight the difference with
gradient-free methods further below.

.. _GRAPE:

GRadient Ascent Pulse Engineering (GRAPE)
-----------------------------------------

The GRAPE method looks at the direct gradient :math:`\partial J_T/\partial
\epsilon_j` with respect to any control parameter :math:`\epsilon_j`.  The
value of the control parameter is then updated in the direction of the
gradient, although in all practical applications, a numerical estimate of the
Hessian :math:`\partial^2 J_T/\partial \epsilon_j \partial \epsilon_{j^\prime}`
should also be included in calculation of the update. The L-BFGS-B quasi-Newton
method :cite:`ByrdSJSC1995,ZhuATMS97` is most commonly used for this purpose.

The control parameter :math:`\epsilon_j` may be the value of a control field in
a particular time interval. When the control field is a discretization of a
time-continuous control, e.g. using the piecewise-constant scheme as in
:numref:`figkrotovscheme`, and for typical functionals, the numerical effort
required for the calculation of the gradient :math:`\partial J_T/\partial
\epsilon_j` is nearly identical to the effort for a single iteration of
Krotov’s method. In both cases, a forward and backward propagation over the
entire time grid is carried out, which determines almost entirely the numerical
effort. This requirement results from the derivative of the complex overlaps
:math:`\tau_k` between the propagated states :math:`\{\ket{\phi_k(T)}\}` and
the target states :math:`\{\ket{\phi_k^{\tgt}}\}`, on which the standard
functionals are based, cf. Eq. :eq:`tauk`. The relevant term in the gradient is
then :cite:`KhanejaJMR05`

.. math::

   \begin{split}
     \frac{\partial \tau_k}{\partial \epsilon_j}
     &= \frac{\partial}{\partial \epsilon_j}
       \big\langle \phi_k^{\tgt} \big\vert
               \Op{U}^{(i)}_{nt-1} \dots \Op{U}^{(i)}_{j} \dots
               \Op{U}^{(i)}_{1} \big\vert \phi_k \big\rangle \\
     &=
       \underbrace{%
           \big\langle \phi_k^{\tgt} \big\vert
             \Op{U}^{(i)}_{nt-1} \dots \Op{U}^{(i)}_{j+1}}_{%
         \bra{\chi^{(i)}_k(t_{j+1})}
        }
         \, \frac{\partial\Op{U}^{(i)}_{j}}{\partial\epsilon_j} \,
        \underbrace{%
          \Op{U}^{(i)}_{j-1} \dots \Op{U}^{(i)}_{1} \big\vert
           \phi_k \big\rangle}_{%
         \ket{\phi^{(i)}_k(t_j)}
        }\,,
     \end{split}

with :math:`\Op{U}^{(i)}_j` the time evolution operator for the time
interval :math:`j`, using the guess controls in iteration :math:`(i)` of
the optimization. We end up with backward-propagated states
:math:`\ket{\chi_k(t_{j+1})}` and forward-propagated states
:math:`\ket{\phi_k(t_j)}`. This is to be compared with the first-order
update equation :eq:`krotov_first_order_update` for Krotov’s method.

In this example of (discretized) time-continuous controls, both GRAPE
and Krotov’s method can generally be used interchangeably. The benefits
of Krotov’s method compared to GRAPE are :cite:`EitanPRA11`:

-  Krotov’s method mathematically guarantees monotonic convergence in
   the continuous limit.

-  Krotov’s method does not require a line search to determine the ideal
   magnitude of the pulse update in the direction of the gradient.

-  The sequential nature of Krotov’s update scheme, where information
   from earlier times enters the update at later times,
   cf. :numref:`figkrotovscheme`, results in
   faster convergence than the concurrent update in
   GRAPE :cite:`MachnesPRA11,JaegerPRA14`. However, this
   advantage disappears as the optimization approaches the
   optimum :cite:`EitanPRA11`.

-  The choice of functional :math:`J_T` in Krotov’s method only enters
   in the boundary condition for the backward-propagated states,
   Eq. :eq:`chi_boundary`, while the update equation stays the same otherwise.
   In contrast, for functionals that do not depend trivially on the overlaps
   :math:`\tau_k = \Braket{\phi_k^\tgt}{\phi_k(T)}`, the evaluation of
   the gradient in GRAPE may deviate significantly from its usual form,
   requiring a problem-specific implementation from scratch.

GRAPE has a significant advantage if the controls are not
time-continuous, but are *physically* piecewise constant (“bang-bang
control”). The calculation of the GRAPE-gradient is unaffected by this,
whereas Krotov’s method can break down when the controls are not
approximately continuous. QuTiP contains an implementation of GRAPE
limited to this use case.

GRAPE (with the second-order derivative estimates by the `L-BFGS-B`_
algorithm) has been shown to converge faster than Krotov’s method when
the optimization is close to the optimum. This is because Krotov’s
method only considers the first-order-gradient with respect to the
control field :cite:`EitanPRA11`, and this derivative
vanishes close to the optimum. This is true even for the Krotov update
with the additional non-linear term, discussed in
the :ref:`SecondOrderKrotov`. There, “second-order” refers to the expansion of
the functional with respect to the states, not to the order of the derivative.

.. _L-BFGS-B: https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html

.. _GrapeInQutip:

GRAPE in QuTiP
--------------

An implementation of GRAPE is included in QuTiP, see the `section on Quantum
Optimal Control in the QuTiP docs`_.  It is used via the
:func:`qutip.control.pulseoptim.optimize_pulse` function.
However, some of the design choices in QuTiP's GRAPE effectively limit
the routine to applications with physically piecewise-constant pulses (where
GRAPE has an advantage over Krotov's method, as discussed in the previous
section).

For discretized time-continuous pulses, the implementation of Krotov's method
in :func:`.optimize_pulses` has the following advantages over
:func:`qutip.control.pulseoptim.optimize_pulse`:

* Krotov's method can optimize for more than one control field at the same time
  (hence the name of the routine :func:`.optimize_pulses` compared to
  :func:`~qutip.control.pulseoptim.optimize_pulse`).
* Krotov's method optimizes a list of :class:`.Objective` instances
  simultaneously. The optimization for multiple simultaneous objectives in
  QuTiP's GRAPE implementation is limited to optimizing a quantum gate. Other
  uses of simultaneous objectives, such as optimizing for robustness, are not
  available.
* Krotov's method can start from an arbitrary set of guess controls. In the
  GRAPE implementation, guess pulses can only be chosen from a specific set of
  options (including "random"). Again, this makes sense for a control field
  that is piecewise constant with relatively few switching points, but is very
  disadvantageous for time-continuous controls.
* Krotov's method has complete flexibility in which propagation method is used
  (via the `propagator` argument to :func:`.optimize_pulses`), while QuTiP's
  GRAPE only allows to choose between fixed number of methods for
  time-propagation. Supplying a problem-specific propagator is not possible.

Thus, QuTiP's GRAPE implementation and the implementation of Krotov's method in
this package complement each other, but will not compare directly.

.. _section on Quantum Optimal Control in the QuTiP docs: http://qutip.org/docs/latest/guide/guide-control.html


Gradient-free optimization
--------------------------

In situations where the controls can be reduced to a relatively small
number of controllable parameters (typically less than 20),
gradient-free optimization becomes feasible. The most straightforward
use case are controls with an analytic shape (maybe due to the
constraints of an experimental setup), with just a few free parameters.
As an example, consider control pulses that are restricted to Gaussian
pulses, so that the only free parameters are the peak amplitude and
pulse width. The control parameters are not required to be parameters of
a time-dependent control, but may also be static parameters in the
Hamiltonian, e.g. the polarization of the laser beams utilized in an
experiment :cite:`HornNJP2018`.

A special case of gradient-free optimization is the Chopped RAndom Basis
(CRAB) method :cite:`DoriaPRL11,CanevaPRA2011`. The essence
of CRAB is in the specific choice of the parametrization in terms of a
low-dimensional *random* basis, as the name implies. Thus, it can be
used when the parametrization is not as “obvious” as in the case of
direct free parameters in the pulse shape discussed above. The
optimization itself is normally performed by Nelder-Mead simplex based
on this parametrization, although any other gradient-free method could
be used as well.

An implementation of CRAB is included in QuTiP, see `QuTiP's documentation of
CRAB`_, and uses the same :func:`qutip.control.pulseoptim.optimize_pulse`
interface as the GRAPE method discussed above (:ref:`GrapeInQutip`) with the
same limitations.

Gradient-free optimization does not require backward propagation, but
only a forward propagation of the initial states and the evaluation of
an arbitrary functional :math:`J_T`. It also does not require the
storage of states. However, the number of iterations can grow extremely
large, especially with an increasing number of control parameters. Thus,
an optimization with a gradient-free method is not necessarily more
efficient overall compared to a gradient-based optimization with much
faster convergence. For only a few parameters, however, it can be highly
efficient.

This makes gradient-free optimization useful for “pre-optimization”,
that is, for finding guess controls that are then further optimized with
a gradient-based method :cite:`GoerzEPJQT2015`. A further
benefit of gradient-free optimization is that it can be applied to *any*
functional, even if :math:`\partial J_T / \partial \bra{\phi_k}` or
:math:`\partial J_T / \partial \epsilon_j` cannot be calculated.

A possible drawback of gradient-free optimization is that is also prone
to get stuck in local optimization minima. To some extent, this can be
mitigated by trying different guess pulses, by
re-parametrization :cite:`RachPRA2015`, or by using some of
the *global* methods available in the NLopt
package :cite:`NLOpt`.

Generally, gradient-free optimization can be easily realized directly in
QuTiP or any other software package for the simulation of quantum
dynamics:

-  Write a function that takes an array of optimization parameters as
   input and returns a figure of merit. This function would, e.g.,
   construct a numerical control pulse from the control parameters,
   simulate the dynamics using :func:`qutip.mesolve.mesolve`, and evaluate a
   figure of merit (like the overlap with a target state).

-  Pass the function to :func:`scipy.optimize.minimize` for gradient-free
   optimization.

The implementation in scipy.optimize.minimize allows to choose between
different optimization methods, with Nelder-Mead simplex being the
default. There exist also more advanced methods such as Subplex_ in
NLopt_ :cite:`NLOpt` that may be worth exploring for
improvements in numerical efficiency, and additional functionality such
as support for non-linear constraints.

.. _Subplex: https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#sbplx-based-on-subplex
.. _NLopt: https://nlopt.readthedocs.io/
.. _QuTiP's documentation of CRAB: http://qutip.org/docs/latest/guide/guide-control.html#the-crab-algorithm


.. _choosing-an-optimization-method:

Choosing an optimization method
-------------------------------

.. _figoctdecisiontree:
.. figure:: oct_decision_tree.svg
   :alt: decision tree.
   :width: 100%

   Decision tree for the choice of an optimization method

Whether to use a gradient-free optimization method, GRAPE, or Krotov’s
method depends on the size of the problem, the requirements on the
control pulse, and the optimization functional. Gradient-free methods
should be used if the number of independent control parameters is
smaller than 20, or the functional is of a form that does not allow to
calculate gradients easily. It is always a good idea to use a
gradient-free method to obtain improved guess pulses for use with a
gradient-based method :cite:`GoerzEPJQT2015`.

GRAPE should be used if the control parameters are discrete, such as on
a coarse-grained time grid, and the derivative of :math:`J_T` with
respect to each control parameter is easily computable. Note that the
implementation provided in QuTiP is limited to state-to-state
transitions and quantum gates, even though the method is generally
applicable to a wider range of objectives.

Krotov’s method should be used if the control is near-continuous, and if
the derivative of :math:`J_T` with respect to the states,
Eq. :eq:`chi_boundary`, can be calculated. When
these conditions are met, Krotov’s method gives excellent convergence.
However, as discussed in the section :ref:`GRAPE`, it is
often observed to slow down when getting close to the minimum of
:math:`J_T`, as the first order derivative vanishes close to the
optimum. For the “best of both worlds”, it can be beneficial to switch
from Krotov’s method to GRAPE with L-BFGS-B in the final stage of the
optimization :cite:`MachnesPRA11`. It has also been proposed
to modify Krotov’s method to include information from the
quasi-Hessian :cite:`EitanPRA11`.

The decision tree in :numref:`figoctdecisiontree` can guide the
choice of an optimization method. The key deciding factor between
gradient-free and gradient-based is the number of control parameters.
For gradient-free optimization, CRAB’s random parametrization is useful
for when there is no obviously better parametrization of the control
(e.g., the control is restricted to an analytic pulse shape and we only
want to optimize the free parameters of that pulse shape). For
gradient-based methods, the decision between GRAPE and Krotov depends
mainly on whether the pulses are approximately time-continuous (up to
discretization), or are of bang-bang type.
