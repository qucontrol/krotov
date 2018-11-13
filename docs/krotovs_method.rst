Krotov’s Method
===============

Functionals
-----------

Krotov's method :cite:`KonnovARC99`, adapted to quantum control,
considers one or more quantum systems, with a set of Hamiltonians :math:`\{\Op{H}_k(\{u_i\})\}`
where each Hamiltonian depends (linearly, for now), on a number of
time-continuous controls. It now seeks it find control fields that optimally
steer a set of states :math:`\{\phi_k\}` in some desired way. To this end, it
minimizes a functional of the form

.. math::

   J[\epsilon^{(i)}(t)]
     = J_T(\{\phi_k^{(i)}(T)\})
         + \int_0^T g_a[\epsilon^{(i)}(t)] \mathrm{d} t
         + \int_0^T g_b[\{\phi^{(i)}_k(t)\}] \mathrm{d} t\,.
   \label{eq:J_krotov}

In addition to the final time functional, e.g. :math:`J_{T,\text{sm}}` or
:math:`J_{T,\text{re}}` defined as

.. math::
   :label: JTsm

    J_{T,\text{sm}} = 1- \frac{1}{N^2} \Abs{\sum_{k=1}^{N} \tau_k}^2
            = 1- \frac{1}{N^2} \sum_{k=1}^{N} \sum_{l=1}^{N} \tau_{l}^* \tau_{k}\,,

.. math::
   :label: JTre

    J_{T,\text{re}} = 1-\frac{1}{N} \Re \left[\, \sum_{k=1}^{N} \tau_k \,\right]
            = 1-\frac{1}{N} \sum_{k=1}^{N} \frac{1}{2} \left( \tau_k + \tau_k^* \right)

*running costs* that depend on
the control field and the propagated states at each point in time are
also included. As before, for a gate optimization, the
:math:`\{\ket{\phi^{(i)}_k(t)}\}` are the basis states
:math:`\{\ket{k}\}` propagated under the pulse :math:`\epsilon^{(i)}(t)`
at the current OCT iteration :math:`i`.

Update Equation
---------------

Krotov’s method uses an auxiliary functional to disentangle the
interdependence of the states and the field, allowing to find an updated
:math:`\epsilon^{(i+1)}(t)` such that
:math:`J[\epsilon^{(i+1)}]  < J[\epsilon^{(i)}]` is guaranteed. The
derivation, see Ref. :cite:`ReichJCP12`, yields the
condition

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
:math:`\ket{\phi_k(0)} = \ket{k}` is written as

.. math::
   :label: fw_eqm

   \frac{\partial}{\partial t} \Ket{\phi_k^{(i+1)}(t)}
     = -\frac{\mathrm{i}}{\hbar} \Op{H}^{(i+1)} \Ket{\phi_k^{(i+1)}(t)}\,,

the co-states :math:`\Ket{\chi_k}` are backward-propagated under the
old pulse as

.. math::
   :label: bw_eqm

   \frac{\partial}{\partial t} \Ket{\chi_k^{(i)}(t)}
     = -\frac{\mathrm{i}}{\hbar} \Op{H}^{\dagger\,(i)} \Ket{\chi_k^{(i)}(t)}
       + \left.\frac{\partial g_b}{\partial \Bra{\phi_k}}\right\vert_{\phi^{(i)}(t)}\,,

with the boundary condition

.. math::
   :label: chi_boundary

   \Ket{\chi_k^{(i)}(T)}
      = - \left.\frac{\partial J_T}{\partial \Bra{\phi_k}}\right\vert_{\phi_k^{(i)}(T)}\,.

In Eq. :eq:`krotov_proto_update`, :math:`\sigma(t)` is a scalar function that must be properly
chosen to ensure monotonic convergence. In many cases, it is sufficient
to set :math:`\sigma(t) \equiv 0`, in particular if the equation of
motion is linear (:math:`\Op{H}` does not depend on
:math:`\ket{\phi_k(t)}`), the functional :math:`J_T` is convex, and no
state-dependent constraints are used (:math:`g_b\equiv 0`). Even for
some types of state-dependent constraints :math:`\sigma(t)` may be set
to zero, specifically for keeping the population in an allowed
subspace :cite:`PalaoPRA2008`. However, a state-dependent
constraint adds an inhomogeneity to the equation of motion for
:math:`\ket{\chi_k(t)}`.

Second order update
-------------------

Where :math:`\sigma(t) \neq 0` is required, it can be determined
numerically as shown in Ref. :cite:`ReichJCP12`. In
chapter [chap:pe], final-time functionals that depend higher than
quadratically on the states are considered, while the equation of motion
remains the linear Schrödinger equation. In this case,

.. math::

   \sigma(t) \equiv -\max\left(\varepsilon_A,2A+\varepsilon_A\right)\,,
     \label{eq:sigma_A}

where :math:`\varepsilon_A` is a small non-negative number that can be
used to enforce strict inequality in the second order optimality
condition. The optimal value for :math:`A` in each iteration can be
determined numerically as :cite:`ReichJCP12`

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

In order to obtain an explicit equation for :math:`\epsilon^{(i+1)}(t)`,
a state-dependent running cost :math:`g_a` must be used, and usually
takes the form

.. math::

   g_a[\epsilon(t)]
     = \frac{\lambda_a}{S(t)} \left(\epsilon(t) - \epsilon^{\text{ref}}(t)\right)^2\,,
     \label{eq:g_a_ref}

with a scaling parameter :math:`\lambda_a` and a shape function
:math:`S(t) \in [0,1]`. When :math:`\epsilon^{\text{ref}}` is set to the optimized
field :math:`\epsilon^{(i)}` from the previous iteration,

.. math::

   g_a[\epsilon^{(i+1)}(t)]
     = \frac{\lambda_a}{S(t)} \left(\Delta\epsilon(t)\right)^2\,,
     \quad
     \Delta\epsilon(t) \equiv \epsilon^{(i+1)}(t) - \epsilon^{(i)}(t)\,,
     \label{eq:g_a_delta}

and for :math:`\sigma(t) \equiv 0`, the explicit first-order Krotov
update equation is obtained :cite:`SklarzPRA2002,PalaoPRA2003`,

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
there, and thus a smooth switch-on and switch-off can be maintained. The
scaling factor :math:`\lambda_a` controls the overall magnitude of the
pulse update. Values that are too large will change
:math:`\epsilon^{(i)}(t)` by only a small amount, causing slow
convergence.

The functional :math:`J_T` enters the first-order update equation only
in the boundary condition for the backward propagated co-state, Eq. :eq:`chi_boundary`.
For the standard functionals defined in Eq. :eq:`JTsm` and Eq. :eq:`JTre`, this evaluates to

.. math::

   \begin{aligned}
     - \left.\frac{\partial J_{T,\text{sm}}}{\partial \Bra{\phi_k}}\right\vert_{\phi_k^{(i)}(T)}
    &= \left( \frac{1}{N^2} \sum_{l=1}^N \tau_l \right) \Op{O} \Ket{k}\,,
    \label{eq:chi_JTsm}
    \\
     - \left.\frac{\partial J_{T,\text{re}}}{\partial \Bra{\phi_k}}\right\vert_{\phi_k^{(i)}(T)}
    &= \frac{1}{2N} \Op{O} \Ket{k}\,.
    \label{eq:chi_JTre}\end{aligned}


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


Discretization
--------------

.. _figkrotovscheme:
.. figure:: krotovscheme.svg
   :alt: Sequential update scheme in Krotov’s method on a time grid.
   :width: 100%

   Sequential update scheme in Krotov’s method on a time grid.

Discretization to a time grid yields the numerical scheme shown in
:numref:`figkrotovscheme`, and resolves the seeming contradiction that the
calculation of :math:`\epsilon^{(i+1)}(t)` requires knowledge of the
states :math:`\ket{\Psi_k^{(i+1)}(t)}` propagated under
:math:`\epsilon^{(i+1)}(t)`. The scheme starts with
:math:`\ket{\chi_k(T)}` obtained from Eq. :eq:`chi_boundary`, which is backward-propagated
under Eq. :eq:`bw_eqm`. All backward-propagated states :math:`\ket{\chi(t)}` must be
stored. The first pulse value is updated according to Eq. :eq:`krotov_first_order_update`, using
:math:`\ket{\chi_k(0)}` and the known initial state
:math:`\ket{\Psi_k(0)} = \ket{k}`. Then, :math:`\ket{\Psi_k(0)}` is
forward-propagated by one time step under Eq. :eq:`fw_eqm` using the updated pulse
value. The updates proceed sequentially, until the final
forward-propagated state :math:`\ket{\Psi_k(T)}` is reached. For
numerical stability, it is useful to define the normalized

.. math:: \ket{\Psi_k^{\text{bw}}(T)} = \frac{1}{\Norm{\chi_k}} \ket{\chi_{k}(T)}

and then later multiply again with :math:`\Norm{\chi_k}` when
calculating the pulse update.


Choice of Step Width
--------------------

The monotonic convergence
of Krotov's method is only guaranteed in the continuous limit; a coarse
time step must be compensated by larger values of the step width :math:`\lambda_a`,
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
     \sum_{k} \Norm{\chi_k}_{\infty} \Norm{\psi_k}_{\infty}
     \Norm{\frac{\partial \Op{H}}{\partial \epsilon}}_{\infty}
     \stackrel{!}{\le}
     \Norm{\epsilon^{(i)}(t)}_{\infty}\,.

Since :math:`S(t) \in [0,1]` and :math:`\ket{\psi_k}` is normalized,
the condition for :math:`\lambda_a` becomes

.. math::

   \lambda_a \ge
     \frac{1}{\max\Abs{\epsilon^{(i)}(t)}}
     \left[ \sum_{k} \Norm{\chi_k}_{\infty} \right]
     \Norm{\frac{\partial \Op{H}}{\partial \epsilon}}_{\infty}\,.

From a practical point of view, the best strategy is to start the
optimization with a comparatively large value of :math:`\lambda_a`, and
after a few iterations lower :math:`\lambda_a` as far as possible
without introducing numerical instabilities. The value of
:math:`\lambda_a` may be adjusted dynamically with the rate of
convergence. Generally, the optimal choice of :math:`\lambda_a` requires
some trial and error.


Rotating Wave Approximation
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
OCT iteration. When the RWA is used, the control-pulses are
complex-valued. In this case, the Krotov update equation is valid for
both the real and the imaginary part independently. That is, in the
update for the real part of the pulse, all derivatives are also taken
with respect to only the real part, and likewise for the imaginary part.


Optimization in Liouville Space
-------------------------------

The control equations have been written in the notation of Hilbert
space. However, they are equally valid for a gate optimization in
Liouville space, by replacing states with density matrices,
:math:`\Op{H}` with :math:`\Liouville`, and inner products with
Hilbert-Schmidt products.

References
----------

.. bibliography:: refs.bib
   :style: unsrt
