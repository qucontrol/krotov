Time Propagation
================

The numerical effort involved in the optimization is almost entirely within the
simulation of the system dynamics. In every iteration, for every objective, the
system must be "propagated" once forwards in time, and once backwards in time.

The implementation of this time propagation must be inside the user-supplied
routine `propagator` that is passed to :func:`.optimize_pulses` and calculates
the propagation over a single time step. In particular,
:func:`qutip.mesolve.mesolve` is not automatically used for simulating any
dynamics within the optimization.  The signature for this function must be as
follows::

    propagator(H, state, dt, c_ops, backwards=False)

It receives four positional arguments:

* `H` is the system Hamiltonian, in a nested-list format similar to that used
  by :func:`qutip.mesolve.mesolve`, e.g. for a Hamiltonian
  $\Op{H} = \Op{H}_0 + c \Op{H}_1$, where $c$ is the value of a control field
  at a particular point in time, `propagator` would receive a list ``[H0, [H1,
  c]]`` where ``H0`` and ``H1`` are :class:`qutip.Qobj` operators.
* `state` is :class:`qutip.Qobj` state that should be propagated
* `dt` is the time step (a float). It is always positive, even for
  ``backwards=True``.
* `c_ops` is a list of collapse (Lindblad) operators, where each list elements
  is of the same form as `H`. The list may be empty for unitary dynamics.

The function also receives one keyword argument, `backwards`. If passed as
`True`, the `propagator` should propagate backwards in time, which usually
means using -`dt` instead of `dt`. In the context of Krotov's method, the
backward propagation uses the conjugate Hamiltonian or Liouvillian. However,
the `propagator` routine does not need to be aware of this fact: it will
receive the appropriate `H` and `c_ops`. Thus, it should not do any complex
conjugation of operators or pulse values internally.

The ``krotov`` package includes some examples for routines that may be used as
`propagators` in :mod:`krotov.propagators`. However, these routines have no
guarantees to be either general or efficient. They are included for use in the
:ref:`krotov-example-notebooks` only.

For "production use", it is recommended to supply a problem-specific
`propagator` that is highly optimized for speed. You might consider the
use of Cython_. This is key to minimize the runtime of the optimization.

.. _Cython: https://cython.org
