Krotov’s Method
===============

*The following overview has been adapted from Ref* :cite:`GoerzPhd2015`.

Functionals
-----------

Krotov's method :cite:`KonnovARC99`, applied to quantum control, considers one
or more quantum systems with a set of Hamiltonians
:math:`\{\Op{H}_k(\{\epsilon_l(t)\})\}` where each Hamiltonian depends on a set
of time-continuous controls :math:`\{\epsilon_l(t)\}`. It seeks to find control
fields that optimally steer a set of initial states :math:`\{\ket{\phi_k}\}` in
some desired way. To this end, in each iteration :math:`(i)`, it minimizes a
functional of the form

.. math::

   \begin{split}
   J[\{\ket{\phi_k^{(i)}(t)}\}, \{\epsilon_l^{(i)}(t)\}]
     &= J_T(\{\ket{\phi_k^{(i)}(T)}\}) \\
     &\qquad
         + \sum_l \int_0^T g_a(\epsilon_l^{(i)}(t)) \dd t
         + \int_0^T g_b(\{\phi^{(i)}_k(t)\}) \dd t\,.
   \end{split}

where :math:`\ket{\phi_k^{(i)}(T)}` are the time-evolved states initial states
:math:`\ket{\phi_k}` under the (guess) controls :math:`\{\epsilon^{(i)}_l(t)\}`
of the :math:`i`'th iteration.

The functional consists of three parts:

* A final time functional $J_T$. This is the "main" part of the functional, and
  we can usually think of $J$ as being an auxiliary functional in the
  optimization of $J_T$.

* A running cost on the control fields, $g_a$. As we will see below, specific
  forms of running costs are required to obtain a closed-form update equation.
  The typical form, and the only one we consider here (and that is realized in
  the :mod:`krotov` package) is

  .. math::

      g_a(\epsilon_l(t))
          = \frac{\lambda_{a, l}}{S_l(t)} \Delta\epsilon_l^2(t)\,.

  We introduce two parameters, the (inverse) Krotov "step width"
  $\lambda_{a,l}$ and the shape function $S_l(t)$ which can be used to
  influence desired properties of the optimized controls. $\Delta\epsilon_l(t)$
  is the update of the control in a single iteration of the optimization
  algorithm. It is best to think of this running cost as a technical
  requirement, and not to assign physical meaning to it. Note that as the
  optimization converges, $\Delta \epsilon_l(t) \rightarrow 0$, so that the
  minimization of $J$ is equivalent to the minimizatin of $J_T$ (for :math:`g_b
  \equiv 0`).

* An optional state-dependent running cost, $g_b$, e.g., to penalize population
  in a subspace. This is rarely used, as there are other methods to achieve the
  same effect, like placing artificially high dissipation on a "forbidden"
  subspace.

The most commonly used functionals (cf. :mod:`krotov.functionals`) optimize for
a set of initial states :math:`\{\ket{\phi_k}\}` to evolve to a set of target
states :math:`\{\ket{\phi_k^\tgt}\}`.  The functionals can then be expressed in
terms of the complex overlaps of the final-time states with the target states
under the given control. Thus,

.. math::
   :label: tauk

     \tau_k = \Braket{\phi_k(T)}{\phi_k^\tgt}

in Hilbert space, or

.. math::

     \tau_k
     = \langle\!\langle \Op{\rho}_k(T) \vert \Op{\rho}^{\tgt} \rangle\!\rangle
     \equiv \tr\left[\Op{\rho}^{\dagger}_k(T) \Op{\rho}_k^{\tgt} \right]

in Liouville space.

The following functionals $J_T$ can be formed from these complex overlaps, taking
into account that any optimization functional $J_T$ must be real. They differ by the way
they treat the phases $\varphi_k$ in the physical optimization goal
:math:`\ket{\phi_k(T)} \overset{!}{=} e^{i\varphi_k}\ket{\phi_k^{\tgt}}`
:cite:`PalaoPRA2003`:

* Optimize for simultaneous state-to-state transitions, with completely arbitrary phases $\varphi_k$,

  .. math::
      :label: JTss

      J_{T,\text{ss}} = 1- \frac{1}{N} \sum_{k=1}^{N} \Abs{\tau_k}^2\,,

  cf. :func:`.J_T_ss`.

* Optimize for simultaneous state-to-state transitions, with an arbitrary *global* phase, i.e.,
  $\varphi_k = \varphi_{\text{global}}$ for all $k$ with arbitrary $\varphi_{\text{global}}$,

  .. math::
      :label: JTsm

      J_{T,\text{sm}} = 1- \frac{1}{N^2} \Abs{\sum_{k=1}^{N} \tau_k}^2
              = 1- \frac{1}{N^2} \sum_{k=1}^{N} \sum_{k'=1}^{N} \tau_{k'}^* \tau_{k}\,,

  cf. :func:`.J_T_sm`.

* Optimize for simultaneous state-to-state transitions, with a global phase of zero, i.e.,
  $\varphi_k = 0$ for all $k$,

  .. math::
      :label: JTre

      J_{T,\text{re}} = 1-\frac{1}{N} \Re \left[\, \sum_{k=1}^{N} \tau_k \,\right]\,,


  cf. :func:`.J_T_re`.


Conditions for monotonic convergence
------------------------------------

Krotov's method is based on a rigorously examination of the conditions for
constructing updated fields :math:`\epsilon_l^{(i+1)}(t)` such that
:math:`J(\{\ket{\phi_k^{(i+1)}(t)}\}, \{\epsilon_l^{(i+1)}\})  \leq
J(\{\ket{\phi_k^{(i)}(t)}\}, \{\epsilon_l^{(i)}\})` is mathematically
guaranteed. The main difficulty is disentangling the
interdependence of the states and the field. Krotov tackles
this by introducing an auxiliary functional :math:`L[\{\ket{\phi_k^{(i)}(t)}\},
\{\epsilon_l^{(i)}(t)\}, \Phi]` that is equivalent to
:math:`J[\{\ket{\phi_k^{(i)}(t)}\}, \{\epsilon_l^{(i)}(t)\}]`, but includes an
arbitrary scalar potential $\Phi$. The freedom in this scalar potential is then
used to formulate a condition for monotonic convergence,

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
       \Bigg\rangle\,,

assuming the equation of motion for the forward propagation of
:math:`\ket{\phi_k}` under the optimized controls to be written as

.. math::
   :label: fw_eqm

   \frac{\partial}{\partial t} \Ket{\phi_k^{(i+1)}(t)}
     = -\frac{\mathrm{i}}{\hbar} \Op{H}^{(i+1)} \Ket{\phi_k^{(i+1)}(t)}\,.

The co-states :math:`\Ket{\chi_k^{(i)}(t)}` are propagated backwards under the
guess controls of iteration (i), i.e., the optimized controls from the previous
iteration, as

.. math::
   :label: bw_eqm

   \frac{\partial}{\partial t} \Ket{\chi_k^{(i)}(t)}
     = -\frac{\mathrm{i}}{\hbar} \Op{H}^{\dagger\,(i)} \Ket{\chi_k^{(i)}(t)}
       + \left.\frac{\partial g_b}{\partial \Bra{\phi_k}}\right\vert_{\phi^{(i)}(t)}\,,

with the boundary condition

.. math::
   :label: chi_boundary

   \Ket{\chi_k^{(i)}(T)}
      = - \left.\frac{\partial J_T}{\partial \Bra{\phi_k}}\right\vert_{\phi^{(i)}(T)}\,.

Note that the backward propagation uses the adjoint Hamiltonian, which becomes
relevant for non-Hermitian Hamiltonians or dissipative dynamics in Liouville
space.  In Hilbert space, and without any state-dependent constraints
(:math:`g_b \equiv 0`), this is still the standard Schrödinger equation running
backwards in time (:math:`\dd t \rightarrow -\dd t`). The equations in
Liouville space follow an analogous structure, with :math:`\Op{H} \rightarrow i
\Liouville`, see :mod:`krotov.mu` for details. A state-dependent constraint
introduces an inhomogeneity. For details on the derivation of the above
equations, see Ref. :cite:`ReichJCP12`.  Here, and in the following, we have
dropped the index :math:`l` of the controls and the associated $\lambda_{a,l}$
and $S_l(t)$; all equations are valid for each individual control.


First order update equation
---------------------------

In order to obtain an explicit equation for :math:`\epsilon^{(i+1)}(t)` --
the optimized pulse in iteration :math:`(i)` -- a running cost
:math:`g_a(\epsilon^{(i+1)}(t))` must be specified. It usually
takes the form

.. math::

   g_a(\epsilon^{(i+1)}(t))
     = \frac{\lambda_a}{S(t)} (\epsilon^{(i+1)}(t) - \epsilon^{\text{ref}}(t))^2\,,

with a scaling parameter :math:`\lambda_a` and a shape function
:math:`S(t) \in [0,1]`. When :math:`\epsilon^{\text{ref}}(t)` is set to the guess
pulse :math:`\epsilon^{(i)}(t)` of the iteration :math:`(i)` (the optimized
pulse from the previous iteration), this yields

.. math::

   g_a(\epsilon^{(i+1)}(t))
     = \frac{\lambda_a}{S(t)} \Delta\epsilon^2(t)\,,
     \quad
     \Delta\epsilon(t) \equiv \epsilon^{(i+1)}(t) - \epsilon^{(i)}(t)\,.

Thus, we obtain the first-order Krotov update equation as :cite:`PalaoPRA2003,SklarzPRA2002`,

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
     \right]\,.

If :math:`S(t) \in [0,1]` is chosen as a function that smoothly goes to
zero at :math:`t=0` and :math:`t=T`, then the update will be suppressed
near the edges of the optimization time interval. Thus, a smooth switch-on
and switch-off can be maintained. The
scaling factor :math:`\lambda_a` controls the overall magnitude of the
pulse update thereby taking the role of an (inverse) "step width".
Values that are too large will change
:math:`\epsilon^{(i)}(t)` by only a small amount in every iteration, causing slow
convergence. Values that are too small will cause sharp spikes in the optimized
control, and numerical instabilities (including a loss of monotonic convergence).

We have assumed that the Hamiltonian is linear in the controls. If this is not
the case, :math:`\epsilon^{(i+1)}(t)` will still show up on the right hand side of
Eq. :eq:`krotov_first_order_update`. In order for
Eq. :eq:`krotov_first_order_update` to remain a valid update equation, we
approximate :math:`\epsilon^{(i+1)}(t) \approx \epsilon^{(i)}(t)` on the right
hand side, that is, :math:`\Abs{\Delta \epsilon(t)} \ll \Abs{\epsilon(t)}`.
This can can be ensured by a sufficiently large value for $\lambda_a$.

The functional :math:`J_T` enters the update equation only implicitly in the
boundary condition for the backward propagated co-state,
Eq. :eq:`chi_boundary`.  For example, the standard functionals defined in
Eq. :eq:`JTsm` and Eq. :eq:`JTre` yield

.. math::

   \begin{aligned}
     - \left.\frac{\partial J_{T,\text{sm}}}{\partial \Bra{\phi_k}}\right\vert_{\phi_k^{(i)}(T)}
    &= \left( \frac{1}{N^2} \sum_{l=1}^N \tau_l \right) \Ket{\phi_k^\tgt}\,,
    \\
     - \left.\frac{\partial J_{T,\text{re}}}{\partial \Bra{\phi_k}}\right\vert_{\phi_k^{(i)}(T)}
    &= \frac{1}{2N} \Op{O} \Ket{\phi_k^\tgt}\,,
    \end{aligned}

cf. :func:`.chis_sm`, :func:`.chis_re`.


Second order update equation
----------------------------

The condition :eq:`krotov_proto_update` and the update
Eq. :eq:`krotov_first_order_update` are based on a first-order expansion of the
auxiliary potential $\Phi$ with respect to the states, see
Ref. :cite:`ReichJCP12` for details. This is sufficient in
most cases, in particular if the equation of
motion is linear (:math:`\Op{H}` does not depend on the states
:math:`\ket{\phi_k(t)}`), the functional :math:`J_T` is convex, and no
state-dependent constraints are used (:math:`g_b\equiv 0`). Even for
some types of state-dependent constraints, the first-order expansion is sufficient,
specifically for keeping the population in an allowed
subspace :cite:`PalaoPRA2008`.

When these conditions are not fulfilled, it is still possible to derive
conditions for monotonic convergence via an expansion of $\Phi$ to second order
in the states, resulting in a second term in Eq. :eq:`krotov_proto_update`,

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

.. math:: \ket{\Delta \phi_k(t)} \equiv \ket{\phi_k^{(i+1)}(t)} - \ket{\phi_k^{(i)}(t)}\,.

In Eq. :eq:`krotov_proto_update2`, :math:`\sigma(t)` is a scalar function that must be properly
chosen to ensure monotonic convergence.

As shown in Ref. :cite:`ReichJCP12`, it is possible to numerically approximate
:math:`\sigma(t)`. In Refs :cite:`WattsPRA2015,GoerzPRA2015`, non-convex
final-time functionals that depend higher than
quadratically on the states are considered, for a standard equation of motion
given by a linear Schrödinger equation. In this case,

.. math::

   \sigma(t) \equiv -\max\left(\varepsilon_A,2A+\varepsilon_A\right)\,,
     \label{eq:sigma_A}

where :math:`\varepsilon_A` is a small non-negative number that can be
used to enforce strict inequality in the second order optimality
condition. The optimal value for :math:`A` in each iteration can be
approximated numerically as :cite:`ReichJCP12`

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
   converging results. Since the second order update requires significantly
   more numerical resources (the calculation of the states
   :math:`\ket{\Delta\phi_k(t)}`), you should always try the optimization with
   the first-order update equation first.


Time discretization
-------------------

.. _figkrotovscheme:
.. figure:: krotovscheme.svg
   :alt: Sequential update scheme in Krotov’s method on a time grid.
   :width: 100%

   Sequential update scheme in Krotov’s method on a time grid.


The derivation of Krotov's method assumes time-continuous control fields. In
this case, it mathematically guarantees monotonic convergence. However, for
practical numerical applications, we have to consider controls on a discrete
time grid with $nt$ points running from :math:`t=0` to :math:`t=T`, with a time
step $\dd t$ . The states are defined on the points of the time grid, while the
controls are assumed to be constant on the intervals of the time grid. See the
notebook `Time Discretization in Quantum Optimal Control`_ for details. This
discretization yields the numerical scheme shown in :numref:`figkrotovscheme`.
The scheme proceeds as follows:

1. Construct the states :math:`\ket{\chi_k(T)}` according to
   Eq. :eq:`chi_boundary`. This may depend on the states forward-propagated
   under the optimized pulse from the previous iteration, that is, the guess
   pulse in the current iteration.

2. Perform a backward-propagation using Eq. :eq:`bw_eqm` as the equation of
   motion, over the entire time grid. The resulting state at each point in the
   time grid must be stored in memory.

3. Starting from the known initial state :math:`\ket{\phi_k(t=0)}`, calculate the
   pulse update for the first time step according to
   Eq. :eq:`krotov_first_order_update`, with $t=\dd t/2$ on the left hand side
   (representing the first *interval* in the time grid, on which the control
   pulse is defined), and $t=0$ on the right-hand side (representing the first
   *point* on the time grid). This approximation of :math:`t \approx t + \dd t
   /2` is what constitutes the "time discretization" mathematically, and what
   resolves the seeming contradiction in the time-continuous
   Eq. :eq:`krotov_first_order_update` that the calculation of
   :math:`\epsilon^{(i+1)}(t)` requires knowledge of the states
   :math:`\ket{\phi_k^{(i+1)}(t)}` propagated under
   :math:`\epsilon^{(i+1)}(t)`.

4. Use the updated control field for the first interval to propagate
   :math:`\ket{\phi_k(t=0)} \rightarrow \ket{\phi_k(t=\dd t)}` for a single
   time step, with Eq. :eq:`fw_eqm` as the equation of motion. The updates then
   proceed sequentially, until the final forward-propagated state
   :math:`\ket{\phi_k(T)}` is reached.

   For numerical stability, it is useful to define the normalized states

   .. math::

      \ket{\phi_k^{\text{bw}}(T)} = \frac{1}{\Norm{\ket{\chi_k}}} \ket{\chi_{k}(T)}

   and use those in the backward propagation, and then later multiply again
   with :math:`\Norm{\ket{\chi_k}}` when calculating the pulse update.


Note that for multiple objectives, the scheme can run in parallel, and each
objectives contributes a term to the update, which are then summed. This is the
sum in :eq:`krotov_first_order_update`. See :mod:`krotov.parallelization` for
details. For a second-order update, the forward propagated states from step 4,
both for the current iteration and the previous iteration, must be stored in
memory over the entire time grid.

.. _Time Discretization in Quantum Optimal Control: https://nbviewer.jupyter.org/gist/goerz/21e46ea7b45c9514e460007de14419bd/Krotov_time_discretization.ipynb#


Choice of λₐ
------------

The monotonic convergence of Krotov's method is only guaranteed in the
continuous limit; a coarse time step must be compensated by larger values of
the inverse step width :math:`\lambda_a`, slowing down convergence. Generally,
choosing :math:`\lambda_a` too small will lead to numerical instabilities and
unphysical features in the optimized pulse. A lower limit for :math:`\lambda_a`
can be determined from the requirement that the change
:math:`\Delta\epsilon(t)` should be at most on the same order of magnitude as
the guess pulse :math:`\epsilon^{(i)}(t)` for that iteration. The
Cauchy-Schwarz inequality applied to the update equation  yields

.. math::

   \Norm{\Delta \epsilon(t)}_{\infty}
     \le
     \frac{\Norm{S(t)}}{\lambda_a}
     \sum_{k} \Norm{\ket{\chi_k (t)}}_{\infty} \Norm{\ket{\phi_k (t)}}_{\infty}
     \Norm{\frac{\partial \Op{H}}{\partial \epsilon}}_{\infty}
     \stackrel{!}{\le}
     \Norm{\epsilon^{(i)}(t)}_{\infty}\,,

where :math:`\norm{\partial \Op{H}/\partial \epsilon}_{\infty}` denotes
the supremum norm of the operator norms of the operator
:math:`\partial \Op{H}/\partial \epsilon` obtained at time $t$.
Since :math:`S(t) \in [0,1]` and :math:`\ket{\phi_k}` is normalized,
the condition for :math:`\lambda_a` becomes

.. math::

   \lambda_a \ge
     \frac{1}{\Norm{\epsilon^{(i)}(t)}_{\infty}}
     \left[ \sum_{k} \Norm{\ket{\chi_k(t)}}_{\infty} \right]
     \Norm{\frac{\partial \Op{H}}{\partial \epsilon}}_{\infty}\,.

From a practical point of view, the best strategy is to start the
optimization with a comparatively large value of :math:`\lambda_a`, and
after a few iterations lower :math:`\lambda_a` as far as possible
without introducing numerical instabilities. The value of
:math:`\lambda_a` may be adjusted dynamically with the rate of
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
constructing :math:`\ket{\chi_k(T)}`, and transforming :math:`\ket{\chi_k(T)}`
back to the rotating frame, before starting the backward-propagation.

When the RWA is used, the control fields are
complex-valued. In this case, the Krotov update equation is valid for
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

The control equations have been written in the notation of Hilbert
space. However, they are equally valid for a gate optimization in
Liouville space, by replacing Hilbert space states with density matrices,
:math:`\Op{H}` with :math:`i \Liouville` (cf. :mod:`krotov.mu`), and inner
products with Hilbert-Schmidt products, :math:`\langle  \cdot \vert \cdot
\rangle \rightarrow \langle\!\langle \cdot  \vert \cdot \rangle\!\rangle`, cf.
e.g. Ref :cite:`GoerzNJP2014`.

See the :ref:`/notebooks/04_example_dissipative_qubit_reset.ipynb` for an example.
