"""Routines that can be passed as `propagator` to :func:`.optimize_pulses`"""
from abc import ABC, abstractmethod

import numpy as np
import qutip
import scipy
from qutip.superoperator import mat2vec, vec2mat
from qutip.cy.spmatfuncs import spmvpy_csr
from qutip.cy.spconvert import dense2D_to_fastcsr_fmode

__all__ = ['expm', 'Propagator', 'DensityMatrixODEPropagator']


def expm(H, state, dt, c_ops, backwards=False, initialize=False):
    """Propagate using matrix exponentiation"""
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
            state (qutip.Qobj): The state to propagation
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

    See :class:`qutip.solver.Options` for arguments.
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

    def __call__(
        self, L, rho, dt, c_ops=None, backwards=False, initialize=False
    ):
        """Evaluation of a single propagation step

        Args:
            L (list): A Liouvillian superoperator in qutip's nested-list
                format, with a scalar value in the place of a time-dependency.
                For example, ``[L0, [L1, u]]`` for a drift Liouvillian ``L0``,
                a control Liouvillian ``H1``, and a scalar value ``u`` that is
                a time-dependent control evaluated for a particular point in
                time. If `initialize` is False, only the control values are
                taken into account; any operators are assumed to be identical
                to the internally cached values of `L` during initialization.
            rho (qutip.Qobj): The density matrix to propagate. The passed value
                is ignored unless `initialize` is given as True. Otherwise, it
                is assumed that `rho` matches the (internally stored) state
                that was the result from the previous propagation step.
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
                propagation. This caches `L` (except for the control values)
                and `rho` internally.
        """
        if initialize:
            self._initialize(L, rho, dt, c_ops, backwards)
        else:
            # only update the control values
            for i in self._control_indices:
                self._L_list[i][1] = L[i][1]
        self._t += dt
        self._r.integrate(self._t)
        return qutip.Qobj(
            dense2D_to_fastcsr_fmode(
                vec2mat(self._r.y), rho.shape[0], rho.shape[1]
            ),
            dims=rho.dims,
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
            initial_vector = mat2vec(rho.full()).ravel('F')
        else:
            raise ValueError("rho must be a density matrix")

        # set up the integrator
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


import pdb

class StateHamODEPropagator(Propagator):
    """Propagator for state under a Hamiltonian 

    See :class:`qutip.solver.Options` for arguments.
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
    ):
        self.method = method
        self.order = order
        self.atol = atol
        self.rtol = rtol
        self.nsteps = nsteps
        self.first_step = first_step
        self.min_step = min_step
        self.max_step = max_step
        self._H_list = None  # cached Hamiltonian data
        self._control_indices = None  # which indices in `H` have a control val
        self._r = None  # the integrator
        self._t = 0.0  # time up to which we've integrated

    def __call__(
        self, H, psi, dt, c_ops=None, backwards=False, initialize=False
    ):
        """Evaluation of a single propagation step

        Args:
            H (list): ...
            psi (qutip.Qobj): ...
            dt (float): The time step over which to propagate
            backwards (bool): ...
            initialize (bool): Whether to (re-)initialize for a new
                propagation. This caches `H` (except for the control values)
                and `psi` internally.
        """
        if initialize:
            self._initialize(H, psi, dt, c_ops, backwards)
        else:
            # only update the control values
            for i in self._control_indices:
                self._H_list[i][1] = H[i][1]
        self._t += dt
        self._r.integrate(self._t)
        return qutip.Qobj(self._r.y, dims=psi.dims)

    def _initialize(self, H, psi, dt, c_ops, backwards):
        H_list = []
        control_indices = []
        if not (c_ops is None or len(c_ops) == 0):
            # in principle, we could convert c_ops to a Lindbladian, here
            raise NotImplementedError("c_ops not implemented")
        for (i, spec) in enumerate(H):
            if isinstance(spec, qutip.Qobj):
                h_op = spec
                h_coeff = 1
            elif isinstance(spec, list) and isinstance(spec[0], qutip.Qobj):
                h_op = spec[0]
                h_coeff = spec[1]
                control_indices.append(i)
            else:
                raise ValueError(
                    "Incorrect specification of time-dependent Hamiltonian"
                )
            if h_op.type == 'oper':
                H_list.append([h_op.data, h_coeff, False])
            else:
                raise ValueError(
                    "Incorrect specification of time-dependent Hamiltonian"
                )
        self._H_list = H_list
        self._control_indices = control_indices
        if psi.type == 'ket':
            initial_vector = psi
        else:
            raise ValueError("psi must be a state")
        def _rhs(t, psi, H_list):
            out = np.zeros(psi.shape[0], dtype=complex)
            H = H_list[0][0]
            H_coeff = (-1)**(int(backwards))*-1j*H_list[0][1]
            spmvpy_csr(H.data, H.indices, H.indptr, psi, H_coeff, out)
            for n in range(1, len(H_list)):
                H = H_list[n][0]
                H_coeff = (-1)**(int(backwards))*-1j*H_list[n][1]
                spmvpy_csr(H.data, H.indices, H.indptr, psi, H_coeff, out)
            return out

        # set up the integrator
        r = scipy.integrate.ode(_rhs)
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
        r.set_f_params(self._H_list)
        self._r = r
        self._t = 0.0

