Krotov’s Method
===============

*The following overview has been adapted from Ref* :cite:`GoerzPhd2015`.

Functionals
-----------

Krotov's method :cite:`KonnovARC99`, adapted to quantum control,
considers one or more quantum systems, with a set of Hamiltonians :math:`\{\Op{H}_k(\{\epsilon_l(t)\})\}`
where each Hamiltonian depends on a set of time-continuous controls
:math:`\epsilon_l(t)`. It seeks to find control fields that optimally
steers a set of initial states :math:`\{\ket{\phi_k}\}` in some desired way. To this end, in
each iteration :math:`i`, it minimizes a functional of the form

.. math::

   J[\epsilon_l^{(i)}(t)]
     = J_T(\{\ket{\phi_k^{(i)}(T)}\})
         + \sum_l \int_0^T g_a[\epsilon_l^{(i)}(t)] \mathrm{d} t
         + \int_0^T g_b[\{\phi^{(i)}_k(t)\}] \mathrm{d} t\,.
   \label{eq:J_krotov}

where :math:`\ket{\phi_k^{(i)}(T)}` are the time-evolved states initial states
:math:`\ket{\phi_k}` under the controls :math:`\{\epsilon_l(t)\}` of the
:math:`i`'th iteration.

The functional consists of three parts:

* A final time functional $J_T$. This is the "main" part of the functional, and
  we can usually think of $J$ as being an auxiliary functional in the
  optimization of $J_T$.

* A running cost on the control fields $g_b$. As we will see below, specific forms of
  running costs are required to obtain a closed-form update equation. The
  typical form, and the only one we consider here (and that is realized in the
  ``krotov`` package) is

  .. math::

      g_a[\epsilon(t)]
          = \frac{\lambda_a}{S(t)} \Abs{\Delta\epsilon(t)}^2\,,

  we introduce two parameters, the "inverse Krotov step width" $\lambda_a$ and the
  shape function $S(t)$ which can be used to influence desired properties of
  the optimized controls. $\Delta\epsilon(t)$ is the update of the control in
  a single iteration of the optimization algorithm. It is best to think of
  this running cost as a technical requirement, and not to assign physical
  meaning to it.

* An optional state-dependent running cost $g_b$, e.g., to penalize population
  in a subspace.

The most commonly used functionals optimize for a set of initial states
:math:`\{\ket{\phi_k}\}` to evolve to a set of target states :math:`\{\ket{\phi_k^\tgt}\}`.
The functionals can then be expressed in terms of the complex overlaps of the
final-time states with the target states under the given control. Thus,

.. math::
     \tau_k = \Braket{\phi_k^\tgt}{\phi_k(T)}

in Hilbert space, or

.. math::

   \label{eq:tau_liouville}
     \tau_k = \tr\left[\Op{\rho}_k^{\tgt\,\dagger} \Op{\rho}_k(T) \right]

in Liouville space.

The following functionals $J_T$ can be formed from these complex overlaps, taking
into account that any optimization functional $J_T$ must be real. They differ by the way
they treat the phases $\varphi_k$ in the physical optimization goal
$\ket{\phi_k(T)} \overset{!}{=} e^{i\varphi_k}\ket{\phi_k^{\tgt}}$ :cite:`PalaoPRA2003`:

* Optimize for simultaneous state-to-state transitions, with completely arbitrary phases $\varphi_k$

.. math::
    :label: JTss

    J_{T,\text{ss}} = 1- \frac{1}{N} \sum_{k=1}^{N} \Abs{\tau_k}^2

* Optimize for simultaneous state-to-state transitions, with an arbitrary *global* phase, i.e.,
  $\varphi_k = \varphi_{\text{global}}$ for all $k$ with arbitrary $\varphi_{\text{global}}$

.. math::
    :label: JTsm

    J_{T,\text{sm}} = 1- \frac{1}{N^2} \Abs{\sum_{k=1}^{N} \tau_k}^2
            = 1- \frac{1}{N^2} \sum_{k=1}^{N} \sum_{k'=1}^{N} \tau_{k'}^* \tau_{k}\,,

* Optimize for simultaneous state-to-state transitions, with a global phase of zero, i.e.,
  $\varphi_k = 0$ for all $k$

.. math::
    :label: JTre

    J_{T,\text{re}} = 1-\frac{1}{N} \Re \left[\, \sum_{k=1}^{N} \tau_k \,\right]
            = 1-\frac{1}{N} \sum_{k=1}^{N} \frac{1}{2} \left( \tau_k + \tau_k^* \right)

Note that in state-to-state optimizations in Liouville space the phases $\varphi_k$ do not appear
at all. This is reflected by the fact that the overlaps $\tau_k$ are always real-valued.
Thus the three functionals above coincide when density matrices are employed.


Conditions for the update equation
----------------------------------

Krotov’s method uses an auxiliary functional to disentangle the
interdependence of the states and the field, allowing to find an updated field
:math:`\epsilon^{(i+1)}(t)` such that
:math:`J[\epsilon^{(i+1)}]  \leq J[\epsilon^{(i)}]` is guaranteed.

Here, and in the following, we drop the index :math:`l` of the controls; all equations
are valid for each individual control.

In Hilbert space, Enforcing the monotonicity condition,
:math:`J[\epsilon^{(i+1)}]  \leq J[\epsilon^{(i)}]`,
see Ref. :cite:`ReichJCP12`, yields the
equation

.. math::
   :label: krotov_proto_update

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
    + \right. \\ & \qquad \left.
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

Assuming the equation of motion for the forward propagation of
:math:`\ket{\phi_k}` to be written as

.. math::
   :label: fw_eqm

   \frac{\partial}{\partial t} \Ket{\phi_k^{(i+1)}(t)}
     = -\frac{\mathrm{i}}{\hbar} \Op{H}^{(i+1)} \Ket{\phi_k^{(i+1)}(t)}\,,

the co-states :math:`\Ket{\chi_k}` are propagated backwards under the
old pulse, i.e., the pulse from the previous iteration, as

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

Note that the backward propagation uses the adjoint Hamiltonian which becomes
relevant for non-Hermitian Hamiltonians or dissipative dynamics. The working equations
in Liouville space follow an analogous structure, see, e.g., Ref. :cite:`ReichJCP12`.

In Eq. :eq:`krotov_proto_update`, :math:`\sigma(t)` is a scalar function that must be properly
chosen to ensure monotonic convergence.

First order update equation
---------------------------

In many cases, it is sufficient
to set :math:`\sigma(t) \equiv 0`, in particular if the equation of
motion is linear (:math:`\Op{H}` does not depend on the states
:math:`\ket{\phi_k(t)}`), the functional :math:`J_T` is convex, and no
state-dependent constraints are used (:math:`g_b\equiv 0`). Even for
some types of state-dependent constraints :math:`\sigma(t)` may be set
to zero, specifically for keeping the population in an allowed
subspace :cite:`PalaoPRA2008`. However, a state-dependent
constraint adds an inhomogeneity to the equation of motion for
:math:`\ket{\chi_k(t)}`.

In order to obtain an explicit equation for :math:`\epsilon^{(i+1)}(t)`,
a state-dependent running cost :math:`g_a` must be used. It usually
takes the form

.. math::

   g_a[\epsilon(t)]
     = \frac{\lambda_a}{S(t)} \left(\epsilon(t) - \epsilon^{\text{ref}}(t)\right)^2\,,

with a scaling parameter :math:`\lambda_a` and a shape function
:math:`S(t) \in [0,1]`. When :math:`\epsilon^{\text{ref}}` is set to the optimized
field :math:`\epsilon^{(i)}` from the previous iteration this yields

.. math::

   g_a[\epsilon^{(i+1)}(t)]
     = \frac{\lambda_a}{S(t)} \left(\Delta\epsilon(t)\right)^2\,,
     \quad
     \Delta\epsilon(t) \equiv \epsilon^{(i+1)}(t) - \epsilon^{(i)}(t)\,,

and for :math:`\sigma(t) \equiv 0`, the first-order Krotov
update equation is re-obtained :cite:`PalaoPRA2003,SklarzPRA2002`,

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
pulse update thereby taking the role of an "inverse Krotov step width".
Values that are too large will change
:math:`\epsilon^{(i)}(t)` by only a small amount in every iteration, causing slow
convergence. Values that are too small will cause sharp spikes in the optimized
control, and numerical instabilities (including a loss of monotonic convergence).

The functional :math:`J_T` enters the first-order update equation only
in the boundary condition for the backward propagated co-state, Eq. :eq:`chi_boundary`.
For the standard functionals defined in Eq. :eq:`JTsm` and Eq. :eq:`JTre`, this evaluates to

.. math::

   \begin{aligned}
     \left.\frac{\partial J_{T,\text{sm}}}{\partial \Bra{\phi_k}}\right\vert_{\phi_k^{(i)}(T)}
    &= - \left( \frac{1}{N^2} \sum_{l=1}^N \tau_l \right) \Ket{\phi_k^\tgt}\,,
    \\
     \left.\frac{\partial J_{T,\text{re}}}{\partial \Bra{\phi_k}}\right\vert_{\phi_k^{(i)}(T)}
    &= - \frac{1}{2N} \Op{O} \Ket{\phi_k^\tgt}\,.
    \end{aligned}

Second order update equation
----------------------------

Where :math:`\sigma(t) \neq 0` is required, it can be approximated
numerically as shown in Ref. :cite:`ReichJCP12`. In
Refs :cite:`WattsPRA2015,GoerzPRA2015`, non-convex final-time functionals
that depend higher than
quadratically on the states are considered, for a standard equation of motion
given by linear Schrödinger equation. In this case,

.. math::

   \sigma(t) \equiv -\max\left(\varepsilon_A,2A+\varepsilon_A\right)\,,
     \label{eq:sigma_A}

where :math:`\varepsilon_A` is a small non-negative number that can be
used to enforce strict inequality in the second order optimality
condition. The optimal value for :math:`A` in each iteration can be
approximated numerically as :cite:`ReichJCP12`

.. math::

   A  =
     \frac{2 \sum_{k=1}^{N} \Re\left[
        \langle \chi_k(T) \vert \Delta\phi_k(T) \rangle
     \right]
           + \Delta J_T}
          {\sum_{k=1}^{N} \Abs{\Delta\phi_k(T)}^2}
     \,,

with

.. math:: \Delta J_T \equiv J_T(\{\phi_k^{(i+1)}(T)\}) -J_T(\{\phi_k^{(i)}(T)\})\,.



.. Non-linear Hamiltonians
   -----------------------

..  If :math:`\Op{H}` depends more than linearly on the field, the
    derivative :math:`\left.\frac{\partial \Op{H}}{\partial \epsilon}\right\vert_{{\scriptsize \begin{matrix}\phi^{(i+1)}(t)\\\epsilon^{(i+1)}(t)\end{matrix}}}`
    yields an explicit dependence on :math:`\epsilon^{(i+1)}(t)` on the
    right hand side of Eq. . In this case, the usual approach is to enforce
    :math:`\epsilon^{(i+1)}(t) \approx \epsilon^{(i)}(t)` with a large value
    of :math:`\lambda_a`. Alternatively, :math:`\Delta\epsilon(t)` may be
    determined in a self-consistent loop. This is especially relevant if
    instead of :math:`\epsilon(t)`, a parametrization :math:`\epsilon(u(t))`
    is used, where :math:`u(t)` is the optimized control field. For example,
    :math:`\epsilon(t) = u^2(t)` is used to ensure that
    :math:`\epsilon(t) > 0`, and

..  .. math::

..     \epsilon(t) = \frac{\epsilon_{\max} - \epsilon_{\min}}{2} \tanh(u(t))
                       + \frac{\epsilon_{\max} + \epsilon_{\min}}{2}

..  keeps :math:`\epsilon(t)` bounded between :math:`\epsilon_{\min}` and
    :math:`\epsilon_{\max}` :cite:`MullerQIP11`.


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
time grid.

Discretization to a time grid yields the numerical scheme shown in
:numref:`figkrotovscheme`, and resolves the seeming contradiction that the
calculation of :math:`\epsilon^{(i+1)}(t)` requires knowledge of the
states :math:`\ket{\Psi_k^{(i+1)}(t)}` propagated under
:math:`\epsilon^{(i+1)}(t)`. The scheme starts with
:math:`\ket{\chi_k(T)}` obtained from Eq. :eq:`chi_boundary`, which is backward-propagated
under Eq. :eq:`bw_eqm`. All backward-propagated states :math:`\ket{\chi_k(t)}` must be
stored. The first pulse value is updated according to Eq. :eq:`krotov_first_order_update`, using
:math:`\ket{\chi_k(0)}` and the known initial state
:math:`\ket{\Psi_k(0)}`. Then, :math:`\ket{\Psi_k(0)}` is
forward-propagated by one time step under Eq. :eq:`fw_eqm` using the updated pulse
value. The updates proceed sequentially, until the final
forward-propagated state :math:`\ket{\Psi_k(T)}` is reached. For
numerical stability, it is useful to define the normalized states

.. math:: \ket{\Psi_k^{\text{bw}}(T)} = \frac{1}{\Norm{\ket{\chi_k}}} \ket{\chi_{k}(T)}

and then later multiply again with :math:`\Norm{\ket{\chi_k}}` when
calculating the pulse update.


Choice of λₐ
------------

The monotonic convergence
of Krotov's method is only guaranteed in the continuous limit; a coarse
time step must be compensated by larger values of the inverse step width :math:`\lambda_a`,
slowing down convergence. Generally, choosing :math:`\lambda_a` too
small will lead to numerical instabilities and unphysical features in
the optimized pulse. A lower limit for :math:`\lambda_a` can be
determined from the requirement that the change
:math:`\Delta\epsilon(t)` should be at most on the same order of
magnitude as the guess pulse :math:`\epsilon^{(i)}(t)` for that
iteration. The Cauchy-Schwarz inequality applied to the update equation 
yields

.. math::

   \Norm{\Delta \epsilon(t)}_{\infty}
     \le
     \frac{\Norm{S(t)}}{\lambda_a}
     \sum_{k} \Norm{\ket{\chi_k (t)}}_{\infty} \Norm{\ket{\psi_k (t)}}_{\infty}
     \Norm{\frac{\partial \Op{H}}{\partial \epsilon}}_{\infty}
     \stackrel{!}{\le}
     \Norm{\epsilon^{(i)}(t)}_{\infty}\,,

where $\Norm{\frac{\partial \Op{H}}{\partial \epsilon}}_{\infty}$ denotes
the supremum norm of the operator norms of the operator
$\frac{\partial \Op{H}}{\partial \epsilon}$ obtained at time $t$.
Since :math:`S(t) \in [0,1]` and :math:`\ket{\psi_k}` is normalized,
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

When using the rotating wave approximation (RWA),
it is important to remember that the target
transformation :math:`\Op{O}` is usually defined in the lab frame, not
in the rotating frame. This is relevant for the construction of
:math:`\ket{\chi_k(T)}`. The easiest approach is to transform the result
of the forward propagation :math:`\ket{\phi_k(T)}` from the rotating
frame to the lab frame, then construct :math:`\ket{\chi_k(T)}` for the
next OCT iteration, and transform :math:`\ket{\chi_k(T)}` back to the
rotating frame, before starting the backward-propagation for the next
OCT iteration. When the RWA is used, the control fields are
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



Optimization in Liouville space
-------------------------------

The control equations have been written in the notation of Hilbert
space. However, they are equally valid for a gate optimization in
Liouville space, by replacing states with density matrices,
:math:`\Op{H}` with :math:`\Liouville`, and inner products with
Hilbert-Schmidt products, see, e.g., Ref. :cite:`ReichJCP12`.

.. .. bibliography:: refs.bib
   :cited:
   :style: unsrt
