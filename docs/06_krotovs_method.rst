Krotov’s Method
===============


Optimization functional
-----------------------

Quantum optimal control methods formalize the problem of finding
"control fields" that achieve some physical objective, using the time
evolution of a quantum system, including possible constraints. The most
direct example is a state-to-state transition, that is, for a known
quantum state at time zero to evolve to a specific target state at final
time :math:`T`, controlling, e.g. a chemical
reaction :cite:`TannorJCP1985`. Another example is the
realization of quantum gates, the building blocks of a quantum computer.
In this case, the states forming a computational basis must transform
according to a unitary transformation :cite:`NielsenChuang`.
The control fields might be the amplitudes of a laser pulse, for the
control of a molecular system, RF fields for nuclear magnetic resonance,
or microwave fields for superconducting circuits. There may be multiple
independent controls involved in the dynamics, such as different color
lasers used in the excitation of a Rydberg atom, or different
polarization components of an electric field.

The quantum control methods build on a rich field of classical control
theory :cite:`BellmanBook,PontryaginBook`. This includes
Krotov’s method :cite:`Krotov.book`, which was originally
formulated to optimize the soft landing of a spacecraft from orbit to
the surface of a planet :cite:`KonnovARC99`, before being
applied to quantum mechanical
problems :cite:`SklarzPRA2002`. Fundamentally, they rely on
the variational principle, that is, the minimization of a functional
:math:`J[\{\ket{\phi_k^{(i)}(t)}\}, \{\epsilon_l^{(i)}(t)\}]` that
includes any required constraints via Lagrange multipliers. The
condition for minimizing :math:`J` is then
:math:`\nabla_{\phi_k, \epsilon_l} J = 0`. In rare cases, the
variational calculus can be solved in closed form, based on Pontryagin’s
maximum principle :cite:`PontryaginBook`. Numerical methods
are required in any other case. These start from an initial guess
control (or set of guess controls, if there are multiple controls), and
calculate an update to these controls that will decrease the value of
the functional. The updated controls then become the guess for the next
iteration of the algorithm, until the value of the functional is
sufficiently small, or convergence is reached.

Mathematically, Krotov’s method, when applied to quantum
systems :cite:`Tannor92,ReichJCP12`, minimizes a functional
of the most general form

.. math::
   :label: functional

     J[\{\ket{\phi_k^{(i)}(t)}\}, \{\epsilon_l^{(i)}(t)\}]
       = J_T(\{\ket{\phi_k^{(i)}(T)}\})
           + \sum_l \int_0^T g_a(\epsilon_l^{(i)}(t)) \dd t
           + \int_0^T g_b(\{\phi^{(i)}_k(t)\}) \dd t\,,

where the :math:`\{\ket{\phi_k^{(i)}(T)}\}` are the time-evolved
initial states :math:`\{\ket{\phi_k}\}` under the (guess) controls
:math:`\{\epsilon^{(i)}_l(t)\}` of the :math:`i`\ ’th iteration. In the
simplest case of a single state-to-state transition, the index :math:`k`
vanishes. For the example of a two-qubit quantum gate,
:math:`\{\ket{\phi_k}\}` would be the logical basis states
:math:`\ket{00}`, :math:`\ket{01}`, :math:`\ket{10}`, and
:math:`\ket{11}`. The sum over :math:`l` vanishes if there is only a
single control. For open system dynamics, the states
:math:`\{\ket{\phi_k}\}` may be density matrices.

The functional consists of three parts:

-  A final-time functional :math:`J_T`. This is the "main" part of the
   functional, and we can usually think of :math:`J` as being an
   auxiliary functional in the optimization of :math:`J_T`.

-  A running cost on the control fields, :math:`g_a`. The most commonly
   used expression (and the only one currently supported by the
   :mod:`krotov` package) is :cite:`PalaoPRA2003`

   .. math::
      :label: g_a

      g_a(\epsilon_l^{(i+1)}(t))
      = \frac{\lambda_{a,l}}{S_l(t)} \Delta\epsilon_l^2(t)\,,
        \quad
      \Delta\epsilon_l(t) \equiv \epsilon_l^{(i+1)}(t) - \epsilon_l^{(i)}(t)\,.

   This introduces two parameters for each control, the (inverse)
   Krotov "step width" :math:`\lambda_{a,l}` and the shape function
   :math:`S_l(t) \in [0, 1]`. :math:`\Delta\epsilon_l(t)` is the update
   of the control in a single iteration of the optimization algorithm.
   As we will see below, :math:`\lambda_{a,l}` determines the overall magnitude
   of :math:`\Delta\epsilon_l(t)`, and :math:`S_l(t)` can be used to ensure
   boundary conditions on :math:`\epsilon^{(i+1)}_l(t)`.

-  An optional state-dependent running cost, :math:`g_b`, can be
   employed, e.g., for penalizing population in a
   subspace :cite:`PalaoPRA2008`. This is rarely used, as
   there are other methods to achieve the same effect, like using a
   non-Hermitian Hamiltonian to remove population from the forbidden
   subspace during the time evolution. Currently, the :mod:`krotov` package
   only supports :math:`g_b \equiv 0`.


The most commonly used final-time functionals (cf. :mod:`krotov.functionals`)
optimize for a set of initial states :math:`\{\ket{\phi_k}\}` to evolve to a
set of target states :math:`\{\ket{\phi_k^\tgt}\}`.  The functionals can then
be expressed in terms of the complex overlaps of the target states with the
final-time states under the given control. Thus,

.. math::
   :label: tauk

     \tau_k = \Braket{\phi_k^\tgt}{\phi_k(T)}

in Hilbert space, or

.. math::

     \tau_k
     = \langle\!\langle \Op{\rho}^{\tgt} \vert \Op{\rho}_k(T) \rangle\!\rangle
     \equiv \tr\left[\Op{\rho}_k^{\tgt\,\dagger} \Op{\rho}_k(T) \right]

in Liouville space.

The following functionals :math:`J_T` can be formed from these complex
overlaps, taking into account that any optimization functional :math:`J_T` must
be real. They differ by the way they treat the phases :math:`\varphi_k` in the
physical optimization goal :math:`\ket{\phi_k(T)} \overset{!}{=}
e^{i\varphi_k}\ket{\phi_k^{\tgt}}` :cite:`PalaoPRA2003`:

* Optimize for simultaneous state-to-state transitions, with completely
  arbitrary phases :math:`\varphi_k`,

  .. math::
      :label: JTss

      J_{T,\text{ss}} = 1- \frac{1}{N} \sum_{k=1}^{N} \Abs{\tau_k}^2\,,

  cf. :func:`.J_T_ss`.

* Optimize for simultaneous state-to-state transitions, with an arbitrary
  *global* phase, i.e., :math:`\varphi_k = \varphi_{\text{global}}` for all
  :math:`k` with arbitrary :math:`\varphi_{\text{global}}`,

  .. math::
      :label: JTsm

      J_{T,\text{sm}} = 1- \frac{1}{N^2} \Abs{\sum_{k=1}^{N} \tau_k}^2
              = 1- \frac{1}{N^2} \sum_{k=1}^{N} \sum_{k'=1}^{N} \tau_{k'}^* \tau_{k}\,,

  cf. :func:`.J_T_sm`.

* Optimize for simultaneous state-to-state transitions, with a global phase of zero, i.e.,
  :math:`\varphi_k = 0` for all :math:`k`,

  .. math::
      :label: JTre

      J_{T,\text{re}} = 1-\frac{1}{N} \Re \left[\, \sum_{k=1}^{N} \tau_k \,\right]\,,


  cf. :func:`.J_T_re`.


Update equation
---------------


Krotov’s method is based on a rigorous examination of the conditions for
calculating the updated fields :math:`\epsilon_l^{(i+1)}(t)` such that
:math:`J(\{\ket{\phi_k^{(i+1)}(t)}\}, \{\epsilon_l^{(i+1)}\}) \leq
J(\{\ket{\phi_k^{(i)}(t)}\}, \{\epsilon_l^{(i)}\})` is true
*by construction* :cite:`Krotov.book,KonnovARC99,PalaoPRA2003,ReichJCP12,SklarzPRA2002`.
It achieves this by adding a vanishing quantity to the functional that
disentangles the implicit dependence of :math:`\{\ket{\phi_k}\}` and
:math:`\{\epsilon_l(t)\}` in the variational calculus. Specifically, the
derivation formulates an auxiliary functional :math:`L[\{\ket{\phi_k^{(i)}(t)}\},
\{\epsilon_l^{(i)}(t)\}, \Phi]` that is equivalent to
:math:`J[\{\ket{\phi_k^{(i)}(t)}\}, \{\epsilon_l^{(i)}(t)\}]`, but includes an
arbitrary scalar potential :math:`\Phi`. The freedom in this scalar potential is then
used to formulate a condition to ensure monotonic convergence,

.. math::
   :label: krotov_proto_update

   \left.\frac{\partial g_a}{\partial \epsilon}\right\vert_{\epsilon^{(i+1)}(t)}
   = 2 \Im
     \sum_{k=1}^{N}
     \Bigg\langle
       \chi_k^{(i)}(t)
     \Bigg\vert
       \Bigg(
          \left.\frac{\partial \Op{H}}{\partial \epsilon}\right\vert_{{\scriptsize \begin{matrix}\phi^{(i+1)}(t)\\\epsilon^{(i+1)}(t)\end{matrix}}}
       \Bigg)
     \Bigg\vert
       \phi_k^{(i+1)}(t)
     \Bigg\rangle\,.

For :math:`g_a` as in Eq. :eq:`g_a`, this condition becomes the Krotov update
equation :cite:`Tannor92,PalaoPRA2003,SklarzPRA2002`,

.. math::
   :label: krotov_first_order_update

   \Delta\epsilon(t)
     =
   \frac{S(t)}{\lambda_a} \Im \left[
     \sum_{k=1}^{N}
     \Bigg\langle
       \chi_k^{(i)}(t)
     \Bigg\vert
       \Bigg(
          \left.\frac{\partial \Op{H}}{\partial \epsilon}\right\vert_{{\scriptsize \begin{matrix}\phi^{(i+1)}(t)\\\epsilon^{(i+1)}(t)\end{matrix}}}
      \Bigg)
     \Bigg\vert
       \phi_k^{(i+1)}(t)
     \Bigg\rangle
   \right]\,,

with the equation of motion for the forward propagation of
:math:`\ket{\phi_k}` under the optimized controls
:math:`\epsilon^{(i+1)}(t)` of the iteration :math:`(i)`,

.. math::
   :label: fw_eqm

   \frac{\partial}{\partial t} \Ket{\phi_k^{(i+1)}(t)}
   = -\frac{\mathrm{i}}{\hbar} \Op{H}^{(i+1)} \Ket{\phi_k^{(i+1)}(t)}\,.

For the moment, we have assumed unitary dynamics; the generalization to
open system dynamics will be discussed later in this section. The
co-states :math:`\ket{\chi_k^{(i)}(t)}` are propagated backwards in time
under the guess controls of iteration :math:`(i)`, i.e., the optimized
controls from the previous iteration, as

.. math::
   :label: bw_eqm

   \frac{\partial}{\partial t} \Ket{\chi_k^{(i)}(t)}
   = -\frac{\mathrm{i}}{\hbar} \Op{H}^{\dagger\,(i)} \Ket{\chi_k^{(i)}(t)}
     + \left.\frac{\partial g_b}{\partial \Bra{\phi_k}}\right\vert_{\phi^{(i)}(t)}\,,

with the boundary condition

.. math::
   :label: chi_boundary

   \Ket{\chi_k^{(i)}(T)}
   = - \left.\frac{\partial J_T}{\partial \Bra{\phi_k}}
     \right\vert_{\phi^{(i)}(T)}\,.

Here, and in the following, we have dropped the index :math:`l` of the
controls and the associated :math:`\lambda_{a,l}` and :math:`S_l(t)`;
all equations are valid for each individual control.

Frequently, the control field :math:`\epsilon(t)` is required to be zero
at :math:`t=0` and :math:`t=T` in order to smoothly switch on and off.
To ensure that the update maintains this behavior,
:math:`S(t) \in [0,1]` is chosen as a function with those same
conditions. A typical example is a :func:`.flattop` function

.. math::

    S(t) = \begin{cases}
      B(t; t_0=0, t_1=2 t_{\text{on}})
        & \text{for} \quad 0 < t < t_{\text{on}} \\
      1 & \text{for} \quad t_{\text{on}} \le t \le T - t_{\text{off}} \\
      B(t; t_0=T-2 t_{\text{off}}, t_1=T)
        & \text{for} \quad T - t_{\text{on}} < t < T\,,
    \end{cases}

with the :func:`.blackman` shape :math:`B(t; t_0, t_1)`, which is similar to a
Gaussian, but exactly zero at :math:`t = t_0, t_1`.

The scaling factor :math:`\lambda_a` controls the overall magnitude of
the pulse update, thereby taking the role of an (inverse) "step size".
Values that are too large will change :math:`\epsilon^{(i)}(t)` by only
a small amount in every iteration, causing slow convergence. Values that
are too small will result in numerical instability, see :ref:`ChoiceOfLambdaA`.

The coupled equations :eq:`krotov_first_order_update`-:eq:`chi_boundary` can be
generalized to open system dynamics by replacing Hilbert space states with
density matrices, :math:`\Op{H}` with :math:`i \Liouville`, and brakets with
Hilbert-Schmidt products, :math:`\langle  \cdot \vert \cdot \rangle \rightarrow
\langle\!\langle \cdot  \vert \cdot \rangle\!\rangle`. In full generality,
:math:`\Op{H}` in Eq. :eq:`krotov_first_order_update` is the operator :math:`H`
on the right-hand side of whatever the equation of motion for the forward
propagation of the states is, written in the form :math:`i \hbar \dot\phi = H
\phi`, cf. Eq. :eq:`fw_eqm`, see :mod:`krotov.mu`. Note also that the backward
propagation Eq. :eq:`bw_eqm` uses the adjoint operator, which is relevant both for a
dissipative Liouvillian :cite:`BartanaJCP93,OhtsukiJCP99,GoerzNJP2014` and a
non-Hermitian Hamiltonian :cite:`MullerQIP11,GoerzQST2018`.


.. _SecondOrderKrotov:

Optimization of non-linear problems or non-convex functionals
-------------------------------------------------------------

The condition :eq:`krotov_proto_update` and the update
Eq. :eq:`krotov_first_order_update` are based on a first-order expansion of the
auxiliary potential :math:`\Phi` with respect to the states.
This first order is sufficient if the equation of motion is linear
(:math:`\Op{H}` does not depend on the states :math:`\ket{\phi_k(t)}`),
the functional :math:`J_T` is convex (all the "standard" functionals for
quantum control are convex), and no state-dependent constraints are used
(:math:`g_b\equiv 0`). When these conditions are not fulfilled, it is
still possible to derive an optimization algorithm with monotonic
convergence via a second term in Eq. :eq:`krotov_proto_update`
:cite:`KonnovARC99,ReichJCP12`,

.. math::
   :label: krotov_proto_update2

   \begin{split}
     \left.\frac{\partial g_a}{\partial \epsilon}\right\vert_{\epsilon^{(i+1)}(t)}
     & =
     2 \Im \left[
       \sum_{k=1}^{N}
       \Bigg\langle
         \chi_k^{(i)}(t)
       \Bigg\vert
         \Bigg(
            \left.\frac{\partial \Op{H}}{\partial \epsilon}\right\vert_{{\scriptsize \begin{matrix}\phi^{(i+1)}(t)\\\epsilon^{(i+1)}(t)\end{matrix}}}
         \Bigg)
       \Bigg\vert
         \phi_k^{(i+1)}(t)
       \Bigg\rangle
    \right. \\ & \qquad \left.
       + \frac{1}{2} \sigma(t)
       \Bigg\langle
         \Delta\phi_k(t)
       \Bigg\vert
         \Bigg(
            \left.\frac{\partial \Op{H}}{\partial \epsilon}\right\vert_{{\scriptsize \begin{matrix}\phi^{(i+1)}(t)\\\epsilon^{(i+1)}(t)\end{matrix}}}
        \Bigg)
       \Bigg\vert
         \phi_k^{(i+1)}(t)
       \Bigg\rangle
     \right]\,,
   \end{split}

with

.. math::

   \ket{\Delta \phi_k(t)} \equiv \ket{\phi_k^{(i+1)}(t)} - \ket{\phi_k^{(i)}(t)}\,.

This second term is the "non-linear" or "second order" contribution.
The corresponding update quation is, assuming Eq. :eq:`g_a`,

.. math::
   :label: krotov_second_order_update

   \begin{split}
   \Delta\epsilon(t)
   & =
   \frac{S(t)}{\lambda_a}  \Im \left[
     \sum_{k=1}^{N}
     \Bigg\langle
       \chi_k^{(i)}(t)
     \Bigg\vert
       \Bigg(
          \left.\frac{\partial \Op{H}}{\partial \epsilon}\right\vert_{{\scriptsize \begin{matrix}\phi^{(i+1)}(t)\\\epsilon^{(i+1)}(t)\end{matrix}}}
       \Bigg)
     \Bigg\vert
       \phi_k^{(i+1)}(t)
     \Bigg\rangle
    \right. \\ & \qquad \qquad \quad \left.
     + \frac{1}{2} \sigma(t)
     \Bigg\langle
       \Delta\phi_k(t)
     \Bigg\vert
       \Bigg(
          \left.\frac{\partial \Op{H}}{\partial \epsilon}\right\vert_{{\scriptsize \begin{matrix}\phi^{(i+1)}(t)\\\epsilon^{(i+1)}(t)\end{matrix}}}
      \Bigg)
     \Bigg\vert
       \phi_k^{(i+1)}(t)
     \Bigg\rangle
   \right]\,.
   \end{split}

The prefactor :math:`\sigma(t)` to the second order update is a scalar function
that must be chosen appropriately to ensure monotonic convergence.

As shown in Ref. :cite:`ReichJCP12`, it is possible to numerically approximate
:math:`\sigma(t)`. In Refs :cite:`WattsPRA2015,GoerzPRA2015`, non-convex
final-time functionals that depend higher than quadratically on the states are
considered, for a standard equation of motion given by a linear Schrödinger
equation. In this case,

.. math::

   \sigma(t) \equiv -\max\left(\varepsilon_A,2A+\varepsilon_A\right)\,,
     \label{eq:sigma_A}

where :math:`\varepsilon_A` is a small non-negative number that can be used to
enforce strict inequality in the second order optimality condition. The optimal
value for :math:`A` in each iteration can be approximated numerically
as :cite:`ReichJCP12`

.. math::

   A  =
     \frac{\sum_{k=1}^{N} 2 \Re\left[
        \langle \chi_k(T) \vert \Delta\phi_k(T) \rangle
     \right]
           + \Delta J_T}
          {\sum_{k=1}^{N} \Abs{\Delta\phi_k(T)}^2}
     \,,

cf. :func:`krotov.second_order.numerical_estimate_A`, with

.. math:: \Delta J_T \equiv J_T(\{\phi_k^{(i+1)}(T)\}) -J_T(\{\phi_k^{(i)}(T)\})\,.


See the :ref:`/notebooks/07_example_PE.ipynb` for an example.

.. Note::

   Even when the second order update equation is mathematically required to
   guarantee monotonic convergence, very often an optimization with the
   first-order update equation :eq:`krotov_first_order_update` will give
   converging results. Since the second order update requires
   more numerical resources (calculation and storage of the states
   :math:`\ket{\Delta\phi_k(t)}`), you should always try the optimization with
   the first-order update equation first.


Time discretization
-------------------

.. _figkrotovscheme:
.. figure:: krotovscheme.svg
   :alt: Sequential update scheme in Krotov’s method on a time grid.
   :width: 100%

   Sequential update scheme in Krotov’s method on a time grid.


The derivation of Krotov’s method assumes time-continuous control fields. Only
in this case, monotonic convergence is mathematically guaranteed. However, for
practical numerical applications, we have to consider controls on a discrete
time grid with :math:`nt` points running from :math:`t=0` to :math:`t=T`, with
a time step :math:`\dd t`. The states are defined on the points of the time
grid, while the controls are assumed to be constant on the intervals of the
time grid.  See the notebook `Time Discretization in Quantum Optimal Control`_
for details.  This discretization yields the numerical scheme shown in
:numref:`figkrotovscheme`.  It proceeds as follows :cite:`PalaoPRA2003`:

1. Construct the states :math:`\ket{\chi^{(i)}_k(T)}` according to
   Eq. :eq:`chi_boundary`. These typically
   depend on the states :math:`\{\ket{\phi^{(i)}_k(T)}\}`
   forward-propagated under the optimized pulse from the previous
   iteration, that is, the guess pulse in the current iteration.

2. Perform a backward-propagation using Eq. :eq:`bw_eqm` as the equation of
   motion over the entire time grid. The resulting state at each point in the
   time grid must be stored in memory.

3. Starting from the known initial states
   :math:`\ket{\phi_k} = \ket{\phi_k(t=0)}`, calculate the pulse update
   for the first time step according to Eq. :eq:`krotov_first_order_update`,
   with :math:`t=\dd t/2` on the left-hand side (representing the first
   *interval* in the time grid, on which the control pulse is defined),
   and :math:`t=0` on the right-hand side (representing the first
   *point* on the time grid). This approximation of
   :math:`t \approx t + \dd t /2` is what constitutes the "time
   discretization" mathematically. It resolves the seeming contradiction
   in the time-continuous Eq. :eq:`krotov_first_order_update`
   that the calculation of :math:`\epsilon^{(i+1)}(t)` requires
   knowledge of the states :math:`\ket{\phi_k^{(i+1)}(t)}` obtained from
   a propagation under :math:`\epsilon^{(i+1)}(t)`.

4. Use the updated field :math:`\epsilon^{(i+1)}(\dd t/2)` for the first
   interval to propagate :math:`\ket{\phi_k(t=0)}` for a single time
   step to :math:`\ket{\phi_k^{(i+1)}(t=\dd t)}`, with
   Eq. :eq:`fw_eqm` as the equation of motion. The
   updates then proceed sequentially, until the final forward-propagated
   state :math:`\ket{\phi^{(i+1)}_k(T)}` is reached.

5. The updated control field becomes the guess control for the next
   iteration of the algorithm, starting again at step 1. The
   optimization continues until the value of the functional :math:`J_T`
   falls below some predefined threshold, or convergence is reached,
   i.e., :math:`\Delta J_T` approaches zero so that no further
   significant improvement of :math:`J_T` is to be expected.


For multiple objectives, the scheme can run in parallel, and each objective
contributes a term to the update. Summation of these terms yields the sum
in Eq. :eq:`krotov_first_order_update`. See :mod:`krotov.parallelization` for
details. For a second-order update, the forward propagated states from step 4,
both for the current iteration and the previous iteration, must be stored in
memory over the entire time grid.

.. _Time Discretization in Quantum Optimal Control: https://nbviewer.jupyter.org/gist/goerz/21e46ea7b45c9514e460007de14419bd/Krotov_time_discretization.ipynb#


.. _ChoiceOfLambdaA:

Choice of λₐ
------------

The monotonic convergence of Krotov's method is only guaranteed in the
continuous limit; a coarse time step must be compensated by larger values of
the inverse step size :math:`\lambda_a`, slowing down convergence. Generally,
choosing :math:`\lambda_a` too small will lead to numerical instabilities and
unphysical features in the optimized pulse. A lower limit for :math:`\lambda_a`
can be determined from the requirement that the change
:math:`\Delta\epsilon(t)` should be at most of the same order of magnitude as
the guess pulse :math:`\epsilon^{(i)}(t)` for that iteration. The
Cauchy-Schwarz inequality applied to the update equation yields

.. math::

   \Norm{\Delta \epsilon(t)}_{\infty}
     \le
     \frac{\Norm{S(t)}}{\lambda_a}
     \sum_{k} \Norm{\ket{\chi_k (t)}}_{\infty} \Norm{\ket{\phi_k (t)}}_{\infty}
     \Norm{\frac{\partial \Op{H}}{\partial \epsilon}}_{\infty}
     \stackrel{!}{\le}
     \Norm{\epsilon^{(i)}(t)}_{\infty}\,,

where :math:`\norm{\partial \Op{H}/\partial \epsilon}_{\infty}` denotes the
supremum norm (with respect to time) of the operator norms of the operators
:math:`\partial \Op{H}/\partial \epsilon` obtained at time :math:`t`.  Since
:math:`S(t) \in [0,1]` and :math:`\ket{\phi_k}` is normalized, the condition
for :math:`\lambda_a` becomes

.. math::

   \lambda_a \ge
     \frac{1}{\Norm{\epsilon^{(i)}(t)}_{\infty}}
     \left[ \sum_{k} \Norm{\ket{\chi_k(t)}}_{\infty} \right]
     \Norm{\frac{\partial \Op{H}}{\partial \epsilon}}_{\infty}\,.

From a practical point of view, the best strategy is to start the
optimization with a comparatively large value of :math:`\lambda_a`, and
after a few iterations lower :math:`\lambda_a` as far as possible
without introducing numerical instabilities. The value of
:math:`\lambda_a` may be adjusted dynamically with respect to the rate of
convergence. Generally, the optimal choice of :math:`\lambda_a` requires
some trial and error.


Rotating wave approximation
---------------------------

When using the rotating wave approximation (RWA), it is important to remember
that the target states are usually defined in the lab frame, not in the
rotating frame. This is relevant for the construction of
:math:`\ket{\chi_k(T)}`. When doing a simple optimization, such as a
state-to-state or a gate optimization, the  easiest approach is to transform
the target states to the rotating frame before calculating
:math:`\ket{\chi_k(T)}`. This is both straightforward and numerically
efficient.

Another solution would be to transform the result of the forward propagation
:math:`\ket{\phi_k(T)}` from the rotating frame to the lab frame, then
constructing :math:`\ket{\chi_k(T)}`, and finally to transform
:math:`\ket{\chi_k(T)}` back to the rotating frame, before starting the
backward propagation.

When the RWA is used the control fields are
complex-valued. In this case the Krotov update equation is valid for
both the real and the imaginary part independently. The most straightforward
implementation of the method is for real controls only, requiring that any
complex control Hamiltonian is rewritten as two independent control
Hamiltonians, one for the real part and one for the imaginary part of the
control field. For example,

.. math::

    \epsilon^*(t) \Op{a} + \epsilon(t) \Op{a}^\dagger
    =  \epsilon_{\text{re}}(t) (\Op{a} + \Op{a}^\dagger) + \epsilon_{\text{im}}(t) (i \Op{a}^\dagger - i \Op{a})

with two independent control fields :math:`\epsilon_{\text{re}}(t)= \Re[\epsilon(t)]` and
:math:`\epsilon_{\text{im}}(t) = \Im[\epsilon(t)]`.

See the :ref:`/notebooks/02_example_lambda_system_rwa_complex_pulse.ipynb` for an
example.

Optimization in Liouville space
-------------------------------

The control equations have been written in the notation of Hilbert space.
However, they are equally valid for a gate optimization in Liouville space, by
replacing Hilbert space states with density matrices, :math:`\Op{H}` with
:math:`i \Liouville` (cf. :mod:`krotov.mu`), and inner products with
Hilbert-Schmidt products, :math:`\langle  \cdot \vert \cdot \rangle \rightarrow
\langle\!\langle \cdot  \vert \cdot \rangle\!\rangle`, cf., e.g.,
Ref. :cite:`GoerzNJP2014`.

See the :ref:`/notebooks/04_example_dissipative_qubit_reset.ipynb` for an
example.
