.. _using-krotov-with-qutip:

Using Krotov with QuTiP
=======================

The ``krotov`` package is designed around `QuTiP`_, a very powerful "Quantum
Toolbox" in Python. This means that all operators and states are expressed as
:class:`qutip.Qobj` quantum objects. The :func:`.optimize_pulses` interface
for Krotov's optimization method is closely linked to the interface of QuTiP's
central :func:`~qutip.mesolve.mesolve` routine for simulating the system
dynamics of a closed or open quantum system. In particular, when setting up an
optimization, the system Hamiltonian should be represented by a nested list.
This is, a Hamiltonian of the form :math:`\Op{H} = \Op{H}_0 + \epsilon(t)
\Op{H}_1` is represented as ``H = [H0, [H1, eps]]`` where ``H0`` and ``H1`` are
:class:`~qutip.Qobj` operators, and ``eps`` is a function with signature
``eps(t, args)``, or an array of pulse values of the length the time grid
(`tlist` parameter in :func:`~qutip.mesolve.mesolve`). The operator can depend
on multiple controls, resulting in ``H = [H0, [H1, eps1], [H2, eps2], ...]``.

The central routine provided by the ``krotov`` package is
:func:`.optimize_pulses`. It takes as input a list of objectives, which are
custom tuples of type :class:`.Objective`. Each objective has an
`initial_state`, which is a :class:`qutip.Qobj` representing a Hilbert space
state or density matrix, an `output_state`, and a Hamiltonian `H` in the
nested-list format described above. For dissipative dynamics, it may also
include a list `c_ops` of collapse (Lindblad) operators, where each Lindblad
operator is a :class:`~qutip.Qobj` operator directly, or, for time-dependent
Lindblad operators, a nested list.

In order to simulate the dynamics of the guess pulse, you can use
:meth:`.Objective.mesolve`.

The optimization routine will automatically extract all controls that it can
find in the objectives (in any of the Hamiltonians, or any of the collapse
operators), and iteratively calculate updates to all controls in order to meet
all `objectives` simultaneously. The result of the optimization will be in the
returned :class:`.Result` object, with a list of the optimized controls in the
`optimized_controls` attribute.
The :attr:`.Result.optimized_objectives` property contains a copy of the
objectives with the `optimized_controls` plugged into the Hamiltonian and
collapse operators. The dynamics under the optimized controls can then again be
simulated through :meth:`.Objective.mesolve`.

While the guess controls that are in the `objectives` on input may be functions
or an array of pulse values on the time grid, the output `optimized_controls`
will always be an array of pulse values.

.. _QuTiP: http://qutip.org
