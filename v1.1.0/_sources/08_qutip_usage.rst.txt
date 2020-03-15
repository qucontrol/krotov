.. _using-krotov-with-qutip:

Using Krotov with QuTiP
=======================

The :mod:`krotov` package is designed around `QuTiP`_, a very powerful "Quantum
Toolbox" in Python. This means that all operators and states are expressed as
:class:`qutip.Qobj` quantum objects. The :func:`.optimize_pulses` interface
for Krotov's optimization method is closely linked to the interface of QuTiP's
central :func:`~qutip.mesolve.mesolve` routine for simulating the system
dynamics of a closed or open quantum system. In particular, when setting up an
optimization, the (time-dependent) system Hamiltonian should be represented by
a nested list.  This is, a Hamiltonian of the form :math:`\Op{H} = \Op{H}_0 +
\epsilon(t) \Op{H}_1` is represented as ``H = [H0, [H1, eps]]`` where ``H0``
and ``H1`` are :class:`~qutip.Qobj` operators, and ``eps`` is a function with
signature ``eps(t, args)``, or an array of control values with the length of the
time grid (`tlist` parameter in :func:`~qutip.mesolve.mesolve`). The operator
can depend on multiple controls, resulting in expressions of the form ``H =
[H0, [H1, eps1], [H2, eps2], ...]``.

The central routine provided by the :mod:`krotov` package is
:func:`.optimize_pulses`. It takes as input a list of objectives, each of which
is an instance of :class:`.Objective`. Each objective has an
:attr:`~.Objective.initial_state`, which is a :class:`qutip.Qobj` representing
a Hilbert space state or density matrix, a :attr:`~.Objective.target` (usually
the target state that the :attr:`~.Objective.initial_state` should evolve into
when the objective is fulfilled), and a Hamiltonian :attr:`~.Objective.H` in
the nested-list format described above. For dissipative dynamics,
:attr:`~.Objective.H` should be a Liouvillian, which can be obtained from the
Hamiltonian and a set of Lindblad operators via
:func:`krotov.objectives.liouvillian`. The Liouvillian again is in nested list
format to express time-dependencies. Alternatively, each objective could also
directly include a list :attr:`~.Objective.c_ops` of collapse (Lindblad)
operators , where each collapse operator is a :class:`~qutip.Qobj` operator.
However, this only makes sense if the time propagation routine takes the
collapse operators into account explicitly, such as in the Monte-Carlo
:func:`~qutip.mcsolve.mcsolve`.  Otherwise, the use of
:attr:`~.Objective.c_ops` is strongly discouraged.

If the control function (``eps`` in the above example) relies on the dict
``args`` for static parameters, those ``args`` can be specified via the
`pulse_options` argument in :func:`.optimize_pulses`. See :ref:`HowtoUseArgs`.

In order to simulate the dynamics of the guess control, you can use
:meth:`.Objective.mesolve`, which delegates to :func:`qutip.mesolve.mesolve`.
There is also a related method :meth:`.Objective.propagate` that uses a
different sampling of the control values, see :mod:`krotov.propagators`.

The optimization routine will automatically extract all controls that it can
find in the objectives, and iteratively calculate updates to all controls in
order to meet all `objectives` simultaneously. The result of the optimization
will be in the returned :class:`.Result` object, with a list of the optimized
controls in :attr:`~.Result.optimized_controls`.
The :obj:`~.Result.optimized_objectives` property contains a copy of the
objectives with the :attr:`~.Result.optimized_controls` plugged into the
Hamiltonian or Liouvillian and/or collapse operators. The dynamics under the
optimized controls can then again be simulated through
:meth:`.Objective.mesolve`.

While the guess controls that are in the `objectives` on input may be
functions, or an array of control values on the time grid, the output
:attr:`~.Result.optimized_controls` will always be an array of control values.

.. _QuTiP: http://qutip.org
