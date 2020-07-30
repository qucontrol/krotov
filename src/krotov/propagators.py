r"""Routines that can be passed as `propagator` to :func:`.optimize_pulses`

The numerical effort involved in the optimization is almost entirely within the
simulation of the system dynamics. In every iteration and for every objective,
the system must be "propagated" once forwards in time and once backwards in
time, see also :mod:`krotov.parallelization`.

The implementation of this time propagation must be inside the user-supplied
routine `propagator` that is passed to :func:`.optimize_pulses` and must
calculate the propagation over a single time step. In particular,
:func:`qutip.mesolve.mesolve` is not automatically used for simulating any
dynamics within the optimization.  The signature for any `propagator`
must be the same as the "reference" :func:`expm` propagator:

    >>> str(inspect.signature(krotov.propagators.expm))
    '(H, state, dt, c_ops=None, backwards=False, initialize=False)'

The arguments are as follows (cf. :class:`Propagator`):

* `H` is the system Hamiltonian or Liouvillian, in a nested-list format similar
  to that used by :func:`qutip.mesolve.mesolve`, e.g., for a Hamiltonian
  $\Op{H} = \Op{H}_0 + c \Op{H}_1$, where $c$ is the value of a control field
  at a particular point in time, `propagator` would receive a list ``[H0, [H1,
  c]]`` where ``H0`` and ``H1`` are :class:`qutip.Qobj` operators.
  The nested-list for `H` used here, with scalar values for the controls, is
  obtained internally from the format used by :func:`~qutip.mesolve.mesolve`,
  with time-dependent controls over the entire time grid, via
  :func:`krotov.conversions.plug_in_pulse_values`.
* `state` is the :class:`qutip.Qobj` state that should be propagated, either a
  Hilbert space state, or a density matrix.
* `dt` is the time step (a float). It is always positive, even for
  ``backwards=True``.
* `c_ops` is None, or a list of collapse (Lindblad) operators, where each list
  element is a :class:`qutip.Qobj` instance (or possibly a nested list, for
  time-dependent Lindblad operators. Note that is generally preferred for `H`
  to be a Liouvillian, for dissipative dynamics.
* `backwards` (:class:`bool`): If passed as `True`, the `propagator` should
  propagate backwards in time. In Hilbert space, this means using -`dt` instead
  of `dt`. In Liouville space, there is no difference between forward and
  backward propagation. In the context of Krotov's method, the backward
  propagation uses the conjugate Hamiltonian, respectively Liouvillian.
  However, the `propagator` routine does not need to be aware of this fact: it
  will receive the appropriate `H` and `c_ops`.
* `initialize` (:class:`bool`): A flag to indicate the beginning of a
  propagation over a time grid. If False in subsequent calls, the `propagator`
  may assume that the input `state` is the result of the previous call to
  `propagator`.

.. warning::
    The routines in this module are provided with no guarantee to be either
    general or efficient. The :func:`expm` propagator is exact to machine
    precision, but generally extremely slow.  For "production use", it is
    recommended to supply a problem-specific `propagator` that is highly
    optimized for speed. You might consider the use of Cython_. This is key to
    minimize the runtime of the optimization.

The `initialize` flag enables "stateful" propagators that cache data between
calls. This can significantly improve numerical efficiency.
:class:`DensityMatrixODEPropagator` is an example for such a propagator. In
general, any stateful `propagator` should be an instance of
:class:`Propagator`.

.. _Cython: https://cython.org
"""
from abc import ABC, abstractmethod

import numpy as np
import qutip
import scipy
import threadpoolctl
from qutip.cy.spconvert import dense2D_to_fastcsr_fmode
from qutip.cy.spmatfuncs import spmvpy_csr
from qutip.superoperator import mat2vec, vec2mat


__all__ = ['expm', 'Propagator', 'DensityMatrixODEPropagator']


def expm(H, state, dt, c_ops=None, backwards=False, initialize=False):
    """Propagate using matrix exponentiation.

    This supports `H` being a Hamiltonian (for a Hilbert space `state`) or a
    Liouvillian (for `state` being a density matrix) in nested-list format.
    Collapse operators `c_ops` are not supported. The propagator is not
    stateful, thus `initialize` is ignored. The matrix exponentiation is
    evaluated in single-threaded mode, to prevent accidental nested
    parallelization.
    """
    if c_ops is None:
        c_ops = []
    if len(c_ops) > 0:
        raise NotImplementedError("Liouville exponentiation not implemented")
    assert isinstance(H, list) and len(H) > 0
    eqm_factor = -1j  # factor in front of H on rhs of the equation of motion
    if isinstance(H[0], list):
        if H[0][1].type == 'super':
            eqm_factor = 1
        if backwards:
            eqm_factor = eqm_factor.conjugate()
        A = (eqm_factor * H[0][1]) * H[0][0]
    else:
        if H[0].type == 'super':
            eqm_factor = 1
        if backwards:
            eqm_factor = eqm_factor.conjugate()
        A = eqm_factor * H[0]
    for part in H[1:]:
        if isinstance(part, list):
            A += (eqm_factor * part[1]) * part[0]
        else:
            A += eqm_factor * part
    ok_types = (state.type == 'oper' and A.type == 'super') or (
        state.type in ['ket', 'bra'] and A.type == 'oper'
    )
    if ok_types:
        with threadpoolctl.threadpool_limits(limits=1):
            return ((A * dt).expm())(state)
    else:
        raise NotImplementedError(
            "Cannot handle argument types A:%s, state:%s"
            % (A.type, state.type)
        )


class Propagator(ABC):
    """Abstract base class for stateful propagators"""

    @abstractmethod
    def __call__(
        self, H, state, dt, c_ops=None, backwards=False, initialize=False
    ):
        """Evaluation of a single propagation step

        Args:
            H (list): A Hamiltonian or Liouvillian in qutip's nested-list
                format, with a scalar value in the place of a time-dependency.
                For example, ``[H0, [H1, u]]`` for a drift Hamiltonian ``H0``,
                a control Hamiltonian ``H1``, and a scalar value ``u`` that is
                a time-dependent control evaluated for a particular point in
                time.
            state (qutip.Qobj): The state to propagate
            dt (float): The time step over which to propagate
            c_ops (list or None): A list of Lindblad operators. Using explicit
                Lindblad operators should be avoided: it is usually more
                efficient to convert them into a Lindbladian, passed as `H`
            backwards (bool): Whether the propagation is forward in time or
                backward in time
            initialize (bool): Whether the propagator should (re-)initialize
                for a new propagation, when the propagator is used to advance
                on a time grid, `initialize` should be passed as True for the
                initial time step (0 to `dt` in a forward propagation, or T to
                T-dt for a backward propagation), and False otherwise.

        Note:
            A propagator may assume the propagation to be "sequential"
            when `initialize` is False. That is, the state to propagate is the
            result of the previous call to the propagator.
        """
        pass


class DensityMatrixODEPropagator(Propagator):
    """Propagator for density matrix evolution under a Lindbladian

    See :class:`qutip.solver.Options` for all arguments except `reentrant`.
    Passing True for the `reentrant` re-initializes the propagator in every
    time step.

    Warning:
        By default, the propagator is not "re-entrant". That is, you cannot use
        more than one instance of :class:`DensityMatrixODEPropagator` in the
        same process at the same time. This limitation is due to
        :class:`scipy.integrate.ode` with the "zvode" integrator not being
        re-entrant. Passing ``reentrant=True`` side-steps this problem by
        re-initializating :class:`scipy.integrate.ode` in every time step. This
        makes it possible to use :class:`DensityMatrixODEPropagator` in the
        optimization of multiple objectives, but creates a significant
        overhead.
    """

    def __init__(
        self,
        method='adams',
        order=12,
        atol=1e-8,
        rtol=1e-6,
        nsteps=1000,
        first_step=0,
        min_step=0,
        max_step=0,
        reentrant=False,
    ):
        self.method = method
        self.order = order
        self.atol = atol
        self.rtol = rtol
        self.nsteps = nsteps
        self.first_step = first_step
        self.min_step = min_step
        self.max_step = max_step
        self._L_list = None  # cached Liouvillian data
        self._control_indices = None  # which indices in `L` have a control val
        self._r = None  # the integrator
        self._t = 0.0  # time up to which we've integrated
        self._y = None  # current vectorized state
        self.reentrant = reentrant

    def __call__(
        self, H, state, dt, c_ops=None, backwards=False, initialize=False
    ):
        """Evaluation of a single propagation step

        Args:
            H (list): A Liouvillian superoperator in qutip's nested-list
                format, with a scalar value in the place of a time-dependency.
                For example, ``[L0, [L1, u]]`` for a drift Liouvillian ``L0``,
                a control Liouvillian ``H1``, and a scalar value ``u`` that is
                a time-dependent control evaluated for a particular point in
                time. If `initialize` is False, only the control values are
                taken into account; any operators are assumed to be identical
                to the internally cached values of `H` during initialization.
            state (qutip.Qobj): The density matrix to propagate. The passed
                value is ignored unless `initialize` is given as True.
                Otherwise, it is assumed that `state` matches the (internally
                stored) state that was the result from the previous propagation
                step.
            dt (float): The time step over which to propagate
            c_ops (list or None): An empty list, or None. Since this propagator
                assumes a full Liouvillian, it cannot be combined with Lindblad
                operators.
            backwards (bool): Whether the propagation is forward in time or
                backward in time. Since the equation of motion for a
                Liouvillian and conjugate Liouvillian is the same, this
                parameter has no effect. Instead, for the backward propagation,
                the conjugate Liouvillian must be passed for `L`.
            initialize (bool): Whether to (re-)initialize for a new
                propagation. This caches `H` (except for the control values)
                and `state` internally.
        """
        # H is really an L, but it's a very bad idea for a subclass not to
        # follow the interface of the parent (including argument names).
        # Internally, however, we'll use L instead of H
        if initialize or self.reentrant:
            self._initialize(H, state, dt, c_ops, backwards)
        else:
            if self.reentrant:
                self._initialize_integrator(self._y)
            # only update the control values
            for i in self._control_indices:
                self._L_list[i][1] = H[i][1]
        self._t += dt
        self._r.integrate(self._t)
        self._y = self._r.y
        return qutip.Qobj(
            dense2D_to_fastcsr_fmode(
                vec2mat(self._y), state.shape[0], state.shape[1]
            ),
            dims=state.dims,
            isherm=True,
        )

    @staticmethod
    def _rhs(t, rho, L_list):
        # _rhs being a staticmethod enables the propagator to be pickled (for
        # parallelization)
        out = np.zeros(rho.shape[0], dtype=complex)
        L = L_list[0][0]
        L_coeff = L_list[0][1]
        spmvpy_csr(L.data, L.indices, L.indptr, rho, L_coeff, out)
        for n in range(1, len(L_list)):
            L = L_list[n][0]
            L_coeff = L_list[n][1]
            spmvpy_csr(L.data, L.indices, L.indptr, rho, L_coeff, out)
        return out

    def _initialize(self, L, rho, dt, c_ops, backwards):
        self._initialize_data(L, rho, dt, c_ops, backwards)
        self._initialize_integrator(self._y)

    def _initialize_data(self, L, rho, dt, c_ops, backwards):
        L_list = []
        control_indices = []
        if not (c_ops is None or len(c_ops) == 0):
            # in principle, we could convert c_ops to a Lindbladian, here
            raise NotImplementedError("c_ops not implemented")
        for (i, spec) in enumerate(L):
            if isinstance(spec, qutip.Qobj):
                l_op = spec
                l_coeff = 1
            elif isinstance(spec, list) and isinstance(spec[0], qutip.Qobj):
                l_op = spec[0]
                l_coeff = spec[1]
                control_indices.append(i)
            else:
                raise ValueError(
                    "Incorrect specification of time-dependent Liouvillian"
                )
            if l_op.type == 'super':
                L_list.append([l_op.data, l_coeff, False])
            else:
                raise ValueError(
                    "Incorrect specification of time-dependent Liouvillian"
                )
        self._L_list = L_list
        self._control_indices = control_indices
        if rho.type == 'oper':
            self._y = mat2vec(rho.full()).ravel('F')  # initial state
        else:
            raise ValueError("rho must be a density matrix")

    def _initialize_integrator(self, initial_vector):
        r = scipy.integrate.ode(self._rhs)
        r.set_integrator(
            'zvode',
            method=self.method,
            order=self.order,
            atol=self.atol,
            rtol=self.rtol,
            nsteps=self.nsteps,
            first_step=self.first_step,
            min_step=self.min_step,
            max_step=self.max_step,
        )
        r.set_initial_value(initial_vector)
        r.set_f_params(self._L_list)
        self._r = r
        self._t = 0.0
