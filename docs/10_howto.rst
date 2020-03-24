How-Tos
=======

.. _HowtoGateOptimization:

How to optimize towards a quantum gate
--------------------------------------

To optimize towards a quantum gate :math:`\Op{O}` in a *closed* quantum system,
set one :class:`.Objective` for each state in the logical basis, with the basis
state :math:`\ket{\phi_k}` as the :attr:`~.Objective.initial_state` and
:math:`\Op{O} \ket{\phi_k}` as the :attr:`~.Objective.target`.

You may use :func:`krotov.gate_objectives() <krotov.objectives.gate_objectives>`
to construct the appropriate list of objectives. See the
:ref:`/notebooks/05_example_transmon_xgate.ipynb` for an example. For more
advanced gate optimizations, also see :ref:`HowtoLIOptimization`,
:ref:`HowtoPEOptimization`, :ref:`HowtoDissipativeOptimization`, and
:ref:`HowtoRobustOptimization`.


How to optimize complex-valued control fields
---------------------------------------------

This implementation of Krotov's method requires real-valued control fields. You
must rewrite your Hamiltonian to contain the real part and the imaginary part
of the field as two independent controls. This is always possible. For example,
for a driven harmonic oscillator in the rotating wave approximation, the
interaction Hamiltonian is given by

.. math::

    \Op{H}_\text{int}
    = \epsilon^*(t) \Op{a} + \epsilon(t) \Op{a}^\dagger
    =  \epsilon_{\text{re}}(t) (\Op{a} + \Op{a}^\dagger) + \epsilon_{\text{im}}(t) (i \Op{a}^\dagger - i \Op{a})\,,

where :math:`\epsilon_{\text{re}}(t)= \Re[\epsilon(t)]` and
:math:`\epsilon_{\text{im}}(t) = \Im[\epsilon(t)]` are considered as two
independent (real-valued) controls.

See the :ref:`/notebooks/02_example_lambda_system_rwa_complex_pulse.ipynb` for an example.


.. _HowtoUseArgs:

How to use `args` in time-dependent control fields
--------------------------------------------------

QuTiP requires that the functions that are used to express time-dependencies
have the signature ``func(t, args)`` where `t` is a scalar value for the time
and `args` is a dict containing values for static parameters, see
`QuTiP's documentation on using the args variable`_. Most of the :mod:`krotov`
package's :ref:`krotov-example-notebooks` use closures or hardcoded values
instead of `args`. For example, in the
:ref:`/notebooks/01_example_simple_state_to_state.ipynb`, the Hamiltonian is
defined as::

   def hamiltonian(omega=1.0, ampl0=0.2):
      """Two-level-system Hamiltonian

      Args:
         omega (float): energy separation of the qubit levels
         ampl0 (float): constant amplitude of the driving field
      """
      H0 = -0.5 * omega * qutip.operators.sigmaz()
      H1 = qutip.operators.sigmax()

      def guess_control(t, args):
         return ampl0 * krotov.shapes.flattop(
               t, t_start=0, t_stop=5, t_rise=0.3, func="blackman"
         )

      return [H0, [H1, guess_control]]

Note how `ampl0` is used in `guess_control` as a closure_ from the surrounding
`hamiltonian` scope, `t_stop` and `t_rise` are hardcoded, and `args` is not
used at all. The function could be rewritten as::

   def guess_control(t, args):
      """Initial control amplitude.

      Args:
         t (float): Time value at which to evaluate the control.
         args (dict): Dictionary containing the value "ampl0" with the
           amplitude of the driving field, "t_stop" with the time at which the
           control shape ends, and "t_rise" for the duration of the
           switch-on/switch-off time.
      """
       return args['ampl0'] * krotov.shapes.flattop(
           t,
           t_start=0,
           t_stop=args['t_stop'],
           t_rise=args['t_rise'],
           func="blackman"
       )

   def hamiltonian(omega=1.0):
       """Two-level-system Hamiltonian

       Args:
           omega (float): energy separation of the qubit levels
       """
       H0 = -0.5 * omega * qutip.operators.sigmaz()
       H1 = qutip.operators.sigmax()

       return [H0, [H1, guess_control]]

   ARGS = dict(ampl0=0.2, t_stop=5, t_rise=0.3)

The `ARGS` must be passed to :func:`.optimize_pulses` via the `pulse_options`
parameter::

   pulse_options = {
      guess_control: dict(lambda_a=5, update_shape=S, args=ARGS)
   }

Both :meth:`.Objective.mesolve` and :meth:`.Objective.propagate` take an
optional `args` dict also.

The `args` in `pulse_options` are used automatically when evaluating the
respective initial guess.  Note that the use of `args` does not extend
to `update_shape`, which is always a function of `t` only.  Any other
parameters in the `update_shape` are best set via :func:`functools.partial`,
see the :ref:`/notebooks/03_example_lambda_system_rwa_non_hermitian.ipynb`.

Compare that example to the
:ref:`/notebooks/02_example_lambda_system_rwa_complex_pulse.ipynb`.
In the latter, the values for the parameters in the control fields and the
Hamiltonian are hardcoded, while in the former, all parameters are centrally
defined in a dict which is passed to the optimization and propagation routines.

.. _QuTiP's documentation on using the args variable: http://qutip.org/docs/latest/guide/dynamics/dynamics-time.html#using-the-args-variable
.. _closure: https://www.learnpython.org/en/Closures



How to stop the optimization when the error crosses some threshold
------------------------------------------------------------------

By default, an optimization stops after a predefined number of iterations
(`iter_stop` parameter in :func:`.optimize_pulses`). However, through the
interplay of the `info_hook` and the `check_convergence` routine  passed to
:func:`.optimize_pulses`, the optimization can be stopped based on the
optimization success or the rate of convergence: The `info_hook` routine should
return the value of the optimization functional or error, which is accessible to
`check_convergence` via the :attr:`.Result.info_vals` attribute, see
:mod:`krotov.convergence` for details.

Generally, you should use the :func:`krotov.info_hooks.print_table` function as
an `info_hook`, which receives a function to evaluate the optimization
functional $J_T$ as a parameter. Then, use
:func:`krotov.convergence.value_below` as a `check_convergence` routine to stop
the optimization when $J_T$ falls below some given threshold.

See the thee :ref:`/notebooks/02_example_lambda_system_rwa_complex_pulse.ipynb` for
an example.


How to exclude a control from the optimization
----------------------------------------------

In order to force the optimization to leave any particular control field
unchanged, set its update shape to :func:`krotov.shapes.zero_shape`
in the `pulse_options` that you pass to :func:`.optimize_pulses`.


How to define a new optimization functional
-------------------------------------------

In order to define a new optimization functional :math:`J_T`:

* Decide on what should go in :attr:`.Objective.target` to best describe the
  *physical* control target. If the control target is reached when the
  :attr:`.Objective.initial_state` evolves to a specific target state under the
  optimal control fields, that target state should be included in
  :attr:`~.Objective.target`.

* Define a function `chi_constructor` that calculates the boundary
  condition for the backward-propagation in Krotov's method,

  .. math::

        \ket{\chi_k(T)} \equiv - \left. \frac{\partial J_T}{\partial \bra{\phi_k(T)}} \right\vert_{\ket{\phi_k(T)}}\,,

  or the equivalent experession in Liouville space. This function should
  calculate the states :math:`\ket{\chi_k}` based  on the forward-propagated
  states :math:`\ket{\phi_k(T)}` and the list of objectives. For convenience,
  when :attr:`~.Objective.target` contains a target state, `chi_constructor`
  will also receive `tau_vals` containing the overlaps :math:`\tau_k =
  \Braket{\phi_k^{\tgt}}{\phi_k(T)}`. See :func:`.chis_re` for an example.

* Optionally, define a function that can be used as an `info_hook`
  in :func:`.optimize_pulses` which returns the value
  :math:`J_T`. This is not required to run an optimization since the
  functional is entirely implicit in `chi_constructor`. However, calculating
  the value of the functional is useful for convergence analysis
  (`check_convergence` in :func:`.optimize_pulses`)

See :mod:`krotov.functionals` for some standard functionals. An example for a
more advanced functional is the :ref:`/notebooks/07_example_PE.ipynb`.


How to penalize population in a forbidden subspace
--------------------------------------------------

In principle, :func:`.optimize_pulses` has a `state_dependent_constraint`.
However, this has some caveats. Most notably, it results in an inhomogeneous
equation of motion, which is currently not implemented.

The recommended "workaround" is to place artificially high dissipation on the
levels in the forbidden subspace. A non-Hermitian Hamiltonian is usually a
good way to realize this. See the
:ref:`/notebooks/03_example_lambda_system_rwa_non_hermitian.ipynb`
for an example.


.. _HowtoLIOptimization:

How to optimize towards a two-qubit gate up to single-qubit corrections
-----------------------------------------------------------------------

On many quantum computing platforms, applying arbitrary single-qubit
gates is easy compared to entangling two-qubit gates. A specific
entangling gate like CNOT is combined with single-qubit gates to form a
universal set of gates. For a given physical system, it can be hard to
know a-priori which entangling gates are easy or even possible to
realize. For example, trapped neutral atoms only allow for the
realization of diagonal two-qubit
gates :cite:`JakschPRL2000,GoerzNJP2014` like CPHASE.
However, the CPHASE gate is "locally equivalent" to CNOT: only
additional single-qubit operations are required to obtain one from the
other. A "local-invariants functional" :cite:`MullerPRA11`
defines an optimization with respect to a such a local equivalence
class, and thus is free to find the specific realization of a two-qubit
gate that is easiest to realize.

Use :func:`krotov.objectives.gate_objectives` with ``local_invariants=True`` in
order to construct a list of objectives suitable for an optimization using the
local-invariant functional :cite:`MullerPRA11`. This optimizes towards a
point in the `Weyl chamber`_.

The |weylchamber package|_ contains the suitable `chi_constructor` routines to
pass to :func:`.optimize_pulses`.

The optimization towards a local equivalence class may require use of the
second-order update equation, see :ref:`SecondOrderUpdate`.



.. _HowtoPEOptimization:

How to optimize towards an arbitrary perfect entangler
------------------------------------------------------

The relevant property of a gate is often its entangling power, and the
requirement for a two-qubit gate in a universal set of gates is that it is a
"perfect entangler". A perfect entangler can produce a maximally entangled
state from a separable input state. Since 85% of all two-qubit gates are
perfect entanglers :cite:`WattsE2013,MuszPRA2013`, a functional that targets an
arbitrary perfect entangler :cite:`WattsPRA2015,GoerzPRA2015` solves the
control problem with the least constraints.

The optimization towards an arbitrary perfect entangler is closely related to
an optimization towards a point in the Weyl chamber
(:ref:`HowtoLIOptimization`): It turns out that
in the geometric representation of the `Weyl chamber`_, all the perfect
entanglers lie within a polyhedron, and we can simply minimize the geometric
distance to the surface of this polyhedron.

Use :func:`krotov.objectives.gate_objectives` with ``gate='PE'`` in
order to construct a list of objectives suitable for an optimization using the
perfect entanglers functional :cite:`WattsPRA2015,GoerzPRA2015`.
This is illustrated in the :ref:`/notebooks/07_example_PE.ipynb`.

Again, the `chi_constructor` is available in the |weylchamber package|_.

Both the optimization towards a local equivalence class and an arbitrary perfect
entangler may require use of the second-order update equation, see
:ref:`SecondOrderUpdate`.

.. |weylchamber package| replace:: ``weylchamber`` package
.. _weylchamber package: https://github.com/qucontrol/weylchamber
.. _Weyl chamber: https://weylchamber.readthedocs.io/en/latest/tutorial.html


.. _HowtoDissipativeOptimization:

How to optimize in a dissipative system
---------------------------------------

To optimize a dissipative system, it is sufficient to set an :class:`.Objective`
with a density matrix for the :attr:`~.Objective.initial_state` and
:attr:`~.Objective.target`, and a Liouvillian in :attr:`.Objective.H`.
See the :ref:`/notebooks/04_example_dissipative_qubit_reset.ipynb` for an
example.

Instead of a Liouvillian, it is also possible to set :attr:`.Objective.H` to
the system Hamiltonian, and :attr:`.Objective.c_ops` to the appropriate
Lindblad operators. However, it is generally much more efficient to use
:func:`krotov.objectives.liouvillian` to convert a time-dependent Hamiltonian
and a list of Lindblad operators into a time-dependent Liouvillian. In either
case, the `propagate` routine passed to :func:`~krotov.optimize.optimize_pulses`
must be aware of and compatible with the convention for the objectives.

Specifically for gate optimization, the routine
:func:`~krotov.objectives.gate_objectives`
can be used to automatically set appropriate objectives for an optimization in
Liouville space. The parameter `liouville_states_set` indicates that the system
dynamics are in Liouville space and sets an appropriate choice of matrices that
track the optimization according to Ref. :cite:`GoerzNJP2014`.
See the :ref:`/notebooks/06_example_3states.ipynb` for an example.

For weak dissipation, it may also be possible to avoid the use of density
matrices altogether, and to instead use a non-Hermitian Hamiltonian. For example, you may
use the effective Hamiltonian from the MCWF method :cite:`PlenioRMP1998`,

.. math::

   \Op{H}_{\text{eff}} = \Op{H} - \frac{i}{2} \sum_k \Op{L}_k^\dagger \Op{L}_k\,,

for the Hermitian Hamiltonian :math:`\Op{H}` and the Lindblad operators
:math:`\Op{L}_k`.  Propagating :math:`\Op{H}_{\text{eff}}` (without quantum
jumps) will lead to a decay in the norm of the state corresponding to how much
dissipation the state is subjected to. Numerically, this will usually increase
the value of the optimization functional (that is, the error). Thus the
optimization can be pushed towards avoiding decoherence, without explicitly
performing the optimization in Liouville space. See the
:ref:`/notebooks/03_example_lambda_system_rwa_non_hermitian.ipynb` for an
example.


.. _HowtoRobustOptimization:

How to optimize for robust pulses
---------------------------------


Control fields can be made robust with respect to variations in the
system by performing an "ensemble
optimization" :cite:`GoerzPRA2014`. The idea is to sample a
representative selection of possible system Hamiltonians, and to
optimize over an average of the entire ensemble. In the functional,
Eq. :eq:`functional`, respectively the update
Eq. :eq:`krotov_first_order_update`,
the index :math:`k` now numbers not only the states, but also different
ensemble Hamiltonians: :math:`\Op{H}(\{\epsilon_l(t)\}) \rightarrow \{\Op{H}_k(\{\epsilon_l(t)\})\}`.

The example considered in Ref. :cite:`GoerzPRA2014` is that
of a CPHASE two-qubit gate on trapped Rydberg atoms. Two classical
fluctuations contribute significantly to the gate error: deviations in
the pulse amplitude (:math:`\Omega = 1` ideally), and fluctuations in
the energy of the Rydberg level (:math:`\Delta_{\text{ryd}} = 0`
ideally). Starting from a set of objectives for the unperturbed system, see
:ref:`HowtoGateOptimization`, :func:`~krotov.objectives.ensemble_objectives`
creates an extended set of objectives that duplicates the original objectives
once for each Hamiltonian from a set perturbed Hamiltonian
:math:`\Op{H}(\Omega \neq 1, \Delta_{\text{ryd}} \neq 0)`.
As shown in Ref. :cite:`GoerzNJP2014`, an optimization over the average of all
these objectives  results in controls that are robust over a wide range of
system perturbations.

A simpler example of an ensemble optimization is
:ref:`/notebooks/08_example_ensemble.ipynb`, which considers a state-to-state
transition in a Lamba-System with a dissipative intermediary state.



.. _HowtoSpectralConstraints:

How to apply spectral constraints
---------------------------------

In principle, Krotov's method can include spectral constraints while
maintaining the guarantee for monotonic convergence :cite:`ReichJMO14` .
However, the calculation of the pulse update with such spectral constraints
requires solving a Fredholm equation of the second kind, which has not yet been
implemented numerically. Thus, the ``krotov`` package does not support this
approach (and no such support is planned).

A "cheap" alternative that usually yields good results is to apply a spectral
filter to the optimized pulses after each iteration. The
:func:`.optimize_pulses` function allows this via the
`modify_params_after_iter` argument.

For example, the following function restricts the spectrum of each pulse to a
given range::

    def apply_spectral_filter(tlist, w0, w1):
       """Spectral filter for real-valued pulses.

       The resulting filter function performs a Fast-Fourier-Transform (FFT) of
       each optimized pulse, and sets spectral components for angular
       frequencies below `w0` or above `w1` to zero. The filtered pulse is then
       the result of the inverse FFT, and multiplying again with the update
       shape for the pulse, to ensure that the filtered pulse still fulfills
       the required boundary conditions.

       Args:
           tlist (numpy.ndarray): Array of time grid values. All pulses must be
               defined on the intervals of this time grid
           w0 (float): The lowest allowed (angular) frequency
           w1 (float): The highest allowed (angular) frequency

       Returns:
           callable: A function that can be passed to
           `modify_params_after_iter` to apply the spectral filter.
       """

        dt = tlist[1] - tlist[0]  # assume equi-distant time grid

        n = len(tlist) - 1  # = len(pulse)
        # remember that pulses are defined on intervals of tlist

        w = np.abs(np.fft.fftfreq(n, d=dt / (2.0 * np.pi)))
        # the normalization factor 2π means that w0 and w1 are angular
        # frequencies, corresponding directly to energies in the Hamiltonian
        # (ħ = 1).

        flt = (w0 <= w) * (w <= w1)
        # flt is the (boolean) filter array, equivalent to an array of values 0
        # and 1

        def _filter(**kwargs):
            # same interface as an `info_hook` function
            pulses = kwargs['optimized_pulses']
            shape_arrays = kwargs['shape_arrays']
            for (pulse, shape) in zip(pulses, shape_arrays):
                spectrum = np.fft.fft(pulse)
                # apply the filter by element-wise multiplication
                spectrum[:] *= flt[:]
                # after the inverse fft, we should also multiply with the
                # update shape function. Otherwise, there is no guarantee that
                # the filtered pulse will be zero at t=0 and t=T (assuming that
                # is what the update shape is supposed to enforce). Also, it is
                # important that we overwrite `pulse` in-place (pulse[:] = ...)
                pulse[:] = np.fft.ifft(spectrum).real * shape

        return _filter

This function is passed to :func:`.optimize_pulses` as e.g.

.. code-block:: python

   modify_params_after_iter=apply_spectral_filter(tlist, 0, 7)

to constrain the spectrum of the pulse to angular frequencies
:math:`\omega \in [0, 7]`.
You may want to explore how such a filter behaves in the example of the
:ref:`/notebooks/05_example_transmon_xgate.ipynb`.

Modifying the optimized pulses "manually" through a
``modify_params_after_iter`` function means that we lose all guarantees of
monotonic convergence. If the optimization with a spectral filter does not
converge, you should increase the value of $\lambda_a$ in the `pulse_options`
that are passed to :func:`.optimize_pulses`. A larger value of $\lambda_a$
results in smaller updates in each iteration. This should also translate into
the filter pulses being closer to the unfiltered pulses, increasing the
probability that the changes due to the filter do not undo the monotonic
convergence. You may also find that the optimization fails if the control
problem physically cannot be solved with controls in the desired spectral
range. Without a good physical intuition, trial and error may be
required.


How to limit the amplitude of the controls
------------------------------------------

Amplitude constraints on the control can be realized indirectly through
parametrization :cite:`MuellerPRA2011`. For example, consider the physical
Hamiltonian :math:`\Op{H} = \Op{H}_0 + \epsilon(t) \Op{H}_1`.

There are several possible parametrizations of :math:`\epsilon(t)`
in terms of an unconstrained function :math:`u(t)`:

* For :math:`\epsilon(t) \ge 0`:

   .. math::

      \epsilon(t) = u^2(t)

* For :math:`0 \le \epsilon(t) < \epsilon_{\max}`:

   .. math::

      \epsilon(t) = \epsilon_{\max} \tanh^2\left(u(t)\right)

* For :math:`\epsilon_{\min} < \epsilon(t) < \epsilon_{\max}`:

   .. math::

      \epsilon(t)
         = \frac{\epsilon_{\max} - \epsilon_{\min}}{2}
              \tanh\left(u(t)\right)
            + \frac{\epsilon_{\max} + \epsilon_{\min}}{2}

Krotov's method can now calculate the update :math:`\Delta u(t)` in each
iteration, and then :math:`\Delta \epsilon(t)` via the above equations.

There is a caveat: In the update equation :eq:`krotov_first_order_update`, we
now have the term

.. math::

   \Bigg(
         \left.\frac{\partial \Op{H}}{\partial u}\right\vert_{{\scriptsize \begin{matrix}\phi^{(i+1)}(t)\\u^{(i+1)}(t)\end{matrix}}}
   \Bigg)
   =
   \Bigg(
         \left.\frac{\partial \epsilon}{\partial u}\frac{\partial \Op{H}}{\partial \epsilon}\right\vert_{{\scriptsize \begin{matrix}\phi^{(i+1)}(t)\\u^{(i+1)}(t)\end{matrix}}}
   \Bigg)

on the right hand side. As the dependendence of :math:`\epsilon(t)` on
:math:`u(t)` is non-linear, we are left with a dependency on the unknown
updated parametrization :math:`u^{(i+1)}(t)`. We resolve this by approximating
:math:`u^{(i+1)}(t) \approx u^{(i)}(t)`, or equivalently :math:`\Delta u(t) \ll
u(t)`, which can be enforced by choosing a sufficiently large value of
:math:`\lambda_a` in the `pulse_options` that are passed to
:func:`.optimize_pulses`.

Currently, the ``krotov`` package does not yet support parametrizations in the
above form, although this is a `planned feature <issue23_>`_.
In the meantime, you could modify the control to fit within the desired
amplitude constaints in the same way as applying spectral constaints, see
:ref:`HowtoSpectralConstraints`.


.. _issue23: https://github.com/qucontrol/krotov/issues/23



How to parallelize the optimization
-----------------------------------

Krotov's method is inherently parallel across different objectives. See
:mod:`krotov.parallelization`, and the
:ref:`/notebooks/05_example_transmon_xgate.ipynb` for an example.

It is exceedingly important to ensure that you do not use any accidental nested
parallelization. The :mod:`numpy` library is often eager to run in a
multi-threaded mode that does not combine well with the process-based
parallelization in :mod:`krotov.parallelization`. See
:ref:`HowtoLimitThreadpool`.


.. _HowtoLimitThreadpool:

How to avoid over-subscribing the CPU when using parallelization
----------------------------------------------------------------

A common caveat of parallelization is that the number of numerically intensive
threads or processes should not be larger than the number of CPUs on the
machine. "Oversubscribing" the CPUs can make a parallelized program run slower
by order of magnitudes compared to a serial program!

One consequence of this realization is that *nested parallelizaton* must be
tightly controlled: If your program used process-based parallelization (and
assuming each process can tax a CPU core at 100%), then you must prevent
multiple threads within each process. Depending on how they were compiled, some
of Python's low-level numerical libraries (:mod:`numpy` in particular) are
eager to run in a multi-threaded mode, and it can be surprisingly difficult to
convince them not to do this. In general, you can
`set environment variables to force low-level numerical code into single-threaded mode`_:

.. code-block:: shell

    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

It may be a good idea to set these variables in your ``.bashrc`` (or the
equivalent for whatever shell you are using), and only change their values when
you specifically want to enable multi-threaded execution. You can sometimes set
these variables inside a Python script or notebook, but you must do so before
importing :mod:`numpy`.

The threadpoolctl_ python package is another alternative of eliminating
unexpected multi-threading. The functions in :mod:`krotov.parallelization` use
this package internally to suppress low-level threads. For example, when using
:func:`krotov.parallelization.parallel_map`, you can expected the execution to
be limited to the given `num_cpus`. Also, :func:`.optimize_pulses` by
defaults limits multi-threading, cf. the `limit_thread_pool` argument. Lastly,
:func:`krotov.propagators.expm` ensures that the matrix exponentiation is
calculated single-threadedly.

Always monitor your processes in a tool like htop_ to watch out for unexpected
CPU usage.

.. _set environment variables to force low-level numerical code into single-threaded mode: https://stackoverflow.com/questions/30791550/limit-number-of-threads-in-numpy/31622299#31622299
.. _threadpoolctl: https://github.com/joblib/threadpoolctl
.. _htop: https://hisham.hm/htop/


.. _HowtoStoreResult:

How to prevent losing an optimization result
--------------------------------------------

Optimizations usually take several hundred to several thousand iterations to
fully converge. Thus, the :func:`.optimize_pulses` routine  may require
significant runtime (often multiple days for large problems). Once an
optimization has completed, you are strongly encouraged to store the result to
disk, using :meth:`.Result.dump`.  You may also consider using
:func:`.dump_result` during the `check_convergence` step to dump the current
state of the optimization to disk at regular intervals. This protects you from
losing work if the optimization is interrupted in any way, like an unexpected
crash.

In order to continue after such a crash, you can restore a :class:`.Result`
object containing the recent state of the optimization using
:meth:`.Result.load` (with the original `objectives` and ``finalize=True`` if
the dump file originates from :func:`.dump_result`). You may then call
:func:`.optimize_pulses` and pass the loaded :class:`.Result` object as
`continue_from`.  The new optimization will start from the most recent
optimized controls as a guess, and continue to count iterations from the
previous result. See :ref:`HowtoContinueOptimization` for further details.


.. _HowtoContinueOptimization:

How to continue from a previous optimization
--------------------------------------------

See :ref:`HowtoStoreResult` for how to continue from an optimization that ended
(crashed) prematurely.  Even when an optimization has completed normally, you
may still want to continue with further iterations -- either because you find
that the original `iter_stop` was insufficient to reach full convergence, or
because you would like to modify some parameters, like the λₐ values for
each control. In this case, you can again call :func:`.optimize_pulses` and
pass the :class:`.Result` object from the previous optimization as
`continue_from`. Note that while you are free to change the `pulse_options`
between the two optimization, the `objectives` must remain the same. The
functional (`chi_constructor`) and the `info_hook` should also remain the same
(otherwise, you may and up with inconsistencies in your :class:`.Result`). The
:class:`.Result` object returned by the second optimization will include all
the data from the first optimization.


How to maximize numerical efficiency
------------------------------------

For systems of non-trivial size, the main numerical effort should be in the
simulation of the system dynamics. Every iteration of Krotov's method requires
a full backward propagation and a full forward propagation of the states associated with each
objective, see :mod:`krotov.propagators`. Therefore, the best numerical
efficiency can be achieved by optimizing the performance of the `propagator`
that is passed to :func:`~krotov.optimize.optimize_pulses`.

One possibility is to implement problem-specific propagators, such as
:class:`krotov.propagators.DensityMatrixODEPropagator`. Going further, you
might consider implementing the propagator with the help of lower-level instructions, e.g.,
by using Cython_.

.. _Cython: https://cython.org


How to deal with the optimization running out of memory
-------------------------------------------------------

Krotov's method requires the storage of at least one set of propagated state
over the entire time grid, for each objective. For the second-order update
equation, up to three sets of stored states per objective may be required. In
particular for larger systems and dynamics in Liouville space, the memory
required for storing these states may be prohibitively expensive.

The :func:`~krotov.optimize.optimize_pulses` accepts a `storage` parameter
to which a constructor for an array-like container can be passed wherein the
propagated states will be stored. It is possible to pass custom out-of-memory
storage objects, such as Dask_ arrays. This may carry a significant penalty in
runtime, however, as states will have to be read from disk, or across the
network.

.. _Dask: http://docs.dask.org/en/latest/


How to avoid the overhead of QuTiP objects
------------------------------------------

If you know what you are doing, it is possible to set up an :class:`.Objective`
without any :class:`qutip.Qobj` instances, using arbitrary low-level objects
instead.  See the :ref:`/notebooks/09_example_numpy.ipynb` for an example.
