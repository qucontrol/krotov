"""Routines for formulating objectives.

Objectives, represented as an :class:`Objective` instance, describe the
*physical* objective of an optimization, e.g. a state-to-state transformation,
or a quantum gate. This is distinct from the *mathematical* formulation of an
optimization functional (:mod:`krotov.functionals`). For the same physical
objective, there are usually several different functionals whose minimization
achieve that objective.
"""
import copy
import itertools
import sys
from collections import defaultdict
from functools import partial

import numpy as np
import qutip
from qutip.solver import Options as QutipSolverOptions
from qutip.solver import Result as QutipSolverResult

from .conversions import (
    _nested_list_shallow_copy,
    control_onto_interval,
    discretize,
    extract_controls,
    extract_controls_mapping,
    plug_in_pulse_values,
)


__all__ = [
    'Objective',
    'gate_objectives',
    'ensemble_objectives',
    'liouvillian',
]


FIX_QUTIP_932 = sys.platform != "linux"
"""Workaround for `QuTiP issue 932`_.

If True, in :meth:`Objective.mesolve`, replace any array controls with an
equivalent function. This results in a signficant slowdown of the propagation,
as it circumvents the use of Cython. Defaults to False on Linux, and True on
any non-Linux system.

.. _QuTiP issue 932: https://github.com/qutip/qutip/issues/932
"""


def _adjoint(op, ignore_errors=False):
    """Calculate adjoint of an object, specifically an operator in QuTiP
    nested-list format.

    Controls are not modified.

    If the adjoint of `op` cannot be calculated, raise a :exc:`ValueError` or
    return `op` unchanged if `ignore_errors` is True.
    """
    if isinstance(op, list):
        adjoint_op = []
        for item in op:
            if isinstance(item, list):
                if len(item) != 2:
                    if ignore_errors:
                        return op
                    else:
                        raise ValueError(
                            "%s is not the in the expected format of the "
                            "two-element list '[operator, control]'" % item
                        )
                adjoint_op.append([_adjoint(item[0]), item[1]])
            else:
                adjoint_op.append(_adjoint(item))
        return adjoint_op
    elif op is None:
        return None
    elif isinstance(op, str):
        return op  # e.g. "PE" target
    else:
        try:
            return op.dag()  # qutip
        except AttributeError:
            try:
                return op.conj().T  # numpy
            except AttributeError:
                try:
                    return op.conjugate().transpose()  # numpy-like
                except AttributeError:
                    if ignore_errors:
                        return op
                    else:
                        raise ValueError("Cannot calculate adjoint of %s" % op)


class Objective:
    """A single objective for optimization with Krotov's method.

    Args:
        initial_state (qutip.Qobj): value for :attr:`initial_state`
        H (qutip.Qobj or list): value for :attr:`H`
        target (qutip.Qobj or None): value for :attr:`target`
        c_ops (list or None): value for :attr:`c_ops`

    Example:

        >>> H0 = - 0.5 * qutip.operators.sigmaz()
        >>> H1 = qutip.operators.sigmax()
        >>> eps = lambda t, args: ampl0
        >>> H = [H0, [H1, eps]]
        >>> krotov.Objective(
        ...     initial_state=qutip.ket('0'), target=qutip.ket('1'), H=H
        ... )
        Objective[|Ψ₀(2)⟩ to |Ψ₁(2)⟩ via [H₀[2,2], [H₁[2,2], u₁(t)]]]

    Raises:
        ValueError: If any arguments have an invalid type or structure. This
            can be surpressed by setting the :attr:`type_checking` class
            attribute to False.

    Note:
        Giving collapse operators via :attr:`c_ops` only makes sense if the
        `propagator` passed to :func:`.optimize_pulses` takes them into account
        explicitly. It is strongly recommended to set :attr:`H` as a Lindblad
        operator instead, see :func:`liouvillian`.

    Attributes:
        H (qutip.Qobj or list): The (time-dependent) Hamiltonian or
            Liouvillian in nested-list format, cf.
            :func:`qutip.mesolve.mesolve`. This includes the control fields.
        initial_state (qutip.Qobj): The initial state, as a Hilbert space
            state, or a density matrix.
        target: An object describing the "target" of the optimization, for the
            dynamics starting from :attr:`initial_state`. Usually, this will be
            the target state (the state into which :attr:`initial_state` should
            evolve). More generally, it can be an arbitrary object meeting the
            conventions of a specific `chi_constructor` function that will be
            passed to :func:`.optimize_pulses`.
        c_ops (list or None): List of collapse operators, cf.
            :func:`~qutip.mesolve.mesolve`, in lieu of :attr:`H` being a
            Liouvillian.

    """

    _counter = defaultdict(int)
    _counter['u{count}(t)'] = 1  # it's nicer to start counting controls from 1
    _count_cache = {}
    _default_attribs = ['initial_state', 'H', 'target', 'c_ops']

    str_use_unicode = True
    """Whether the string representation of an :class:`Objective` may use
    unicode symbols, cf. :meth:`summarize` (class attribute)."""

    type_checking = True
    """By default, instantiating :class:`Objective` with invalid types raises a
    :exc:`ValueError`. Setting this to False disables type checks in the
    initializer, allowing certain advanced use cases such as using plain numpy
    objects instead of QuTiP objects (class attribute)."""

    def __init__(self, *, initial_state, H, target, c_ops=None):
        if c_ops is None:
            c_ops = []
        if self.type_checking:
            if not isinstance(H, (qutip.Qobj, list)):
                raise ValueError(
                    "Invalid H, must be a Qobj, or a nested list, not %s"
                    % H.__class__.__name__
                )
        self.H = H
        if self.type_checking:
            if not isinstance(initial_state, qutip.Qobj):
                raise ValueError(
                    "Invalid initial_state: must be Qobj, not %s"
                    % initial_state.__class__.__name__
                )
        self.initial_state = initial_state
        self.target = target
        if self.type_checking:
            if not isinstance(c_ops, list):
                raise ValueError(
                    "Invalid c_ops: must be a list, not %s"
                    % c_ops.__class__.__name__
                )
        self.c_ops = c_ops

    def __copy__(self):
        # When we use copy.copy(objective), we want a
        # semi-deep copy where nested lists in the Hamiltonian and the c_ops
        # are re-created (copy by value), but non-list elements are copied by
        # reference.
        new_objective = Objective(
            H=_nested_list_shallow_copy(self.H),
            initial_state=self.initial_state,
            target=self.target,
            c_ops=[_nested_list_shallow_copy(c) for c in self.c_ops],
        )
        # restore custom attributes (like 'weight')
        for attr in self.__dict__:
            if attr not in self._default_attribs:
                setattr(new_objective, attr, getattr(self, attr))
        return new_objective

    def __deepcopy__(self, memo):
        new_objective = Objective(
            H=copy.deepcopy(self.H, memo),
            initial_state=copy.deepcopy(self.initial_state, memo),
            target=copy.deepcopy(self.target, memo),
            c_ops=[copy.deepcopy(c, memo) for c in self.c_ops],
        )
        # copy custom attributes (like 'weight')
        for attr in self.__dict__:
            if attr not in self._default_attribs:
                setattr(
                    new_objective,
                    attr,
                    copy.deepcopy(getattr(self, attr), memo),
                )
        return new_objective

    def __eq__(self, other):
        if other.__class__ is self.__class__:
            if self.__dict__.keys() != other.__dict__.keys():
                return False
            for attr in self.__dict__:
                a = getattr(self, attr)
                b = getattr(other, attr)
                if not _recursive_eq(a, b):
                    return False
            return True
        else:
            return NotImplemented

    def __ne__(self, other):  # pragma: nocover
        result = self.__eq__(other)
        if result is NotImplemented:
            return NotImplemented
        else:
            return not result

    def adjoint(self):
        """The :class:`Objective` containing the adjoint of all components.

        This does not affect the controls in :attr:`H`: these are
        assumed to be real-valued. Also, :attr:`.Objective.target` will be left
        unchanged if its adjoint cannot be calculated (if it is not a target
        state).
        """
        adjoint_objective = Objective(
            H=_adjoint(self.H),
            initial_state=_adjoint(self.initial_state),
            target=_adjoint(self.target, ignore_errors=True),
            c_ops=[_adjoint(op) for op in self.c_ops],
        )
        # copy custom attributes (like 'weight')
        for attr in self.__dict__:
            if attr not in self._default_attribs:
                setattr(adjoint_objective, attr, getattr(self, attr))
        return adjoint_objective

    def mesolve(
        self,
        tlist,
        *,
        rho0=None,
        H=None,
        c_ops=None,
        e_ops=None,
        args=None,
        **kwargs
    ):
        """Run :func:`qutip.mesolve.mesolve` on the system of the objective.

        Solve the dynamics for the :attr:`H` and :attr:`c_ops` of the
        objective, starting from the objective's :attr:`initial_state`, by
        delegating to :func:`qutip.mesolve.mesolve`. Both the initial state and
        the dynamical generator for the propagation can be overridden by
        passing `rho0` and `H`/`c_ops`.

        Args:
            tlist (numpy.ndarray): array of time grid points on which the
                states are defined
            rho0 (qutip.Qobj or None): The initial state for the propagation.
                If None, the :attr:`initial_state` attribute is used.
            H (qutip.Qobj or None): The dynamical generator (Hamiltonian or
                Liouvillian) for the propagation. If None, the :attr:`H`
                attribute is used.
            c_ops (list or None): List of collapse (Lindblad) operators. If
                None, the :attr:`c_ops` attribute is used.
            e_ops (list or None): A list of operators whose expectation values
                to calculate, for every point in `tlist`. See
                :func:`qutip.mesolve.mesolve`.
            args (dict or None): dictionary of parameters for time-dependent
                Hamiltonians and collapse operators
            **kwargs: All further arguments will be passed to
                :func:`qutip.mesolve.mesolve`.

        Returns:
            qutip.solver.Result: Result of the propagation, see
            :func:`qutip.mesolve.mesolve` for details.
        """
        if rho0 is None:
            rho0 = self.initial_state
        if e_ops is None:
            e_ops = []
        if H is None:
            H = self.H
        if c_ops is None:
            c_ops = self.c_ops
        if args is None:
            args = {}
        if FIX_QUTIP_932:  # pragma: no cover
            controls = extract_controls([self])
            pulses_mapping = extract_controls_mapping([self], controls)
            mapping = pulses_mapping[0]  # "first objective" (dummy structure)
            H = _plug_in_array_controls_as_func(H, controls, mapping[0], tlist)
            c_ops = [
                _plug_in_array_controls_as_func(
                    c_op, controls, mapping[ic + 1], tlist
                )
                for (ic, c_op) in enumerate(self.c_ops)
            ]

        # local instantations for `options` is to work around
        # https://github.com/qutip/qutip/issues/1061
        options = kwargs.pop('options', QutipSolverOptions())

        return qutip.mesolve(
            H=H,
            rho0=rho0,
            tlist=tlist,
            c_ops=c_ops,
            e_ops=e_ops,
            args=args,
            options=options,
            **kwargs
        )

    def propagate(
        self,
        tlist,
        *,
        propagator,
        rho0=None,
        H=None,
        c_ops=None,
        e_ops=None,
        args=None,
        expect=qutip.expect
    ):
        """Propagate the system of the objective over the entire time grid.

        Solve the dynamics for the `H` and `c_ops` of the objective. If `rho0`
        is not given, the `initial_state` will be propagated. This is similar
        to the :meth:`mesolve` method, but instead of using
        :func:`qutip.mesolve.mesolve`, the `propagate` function is used to go
        between points on the time grid. This function is the same as what is
        passed to :func:`.optimize_pulses`. The crucial difference between this
        and :meth:`mesolve` is in the time discretization convention. While
        :meth:`mesolve` uses piecewise-constant controls centered around the
        values in `tlist` (the control field switches in the middle between two
        points in `tlist`), :meth:`propagate` uses piecewise-constant controls
        on the intervals of `tlist` (the control field switches on the points
        in `tlist`). The function `expect` is used to calculate expecation
        values; it receives two parameters, an operator from `e_ops` and a
        state, and must returnd the expectation value of the operator.

        Comparing the result of :meth:`mesolve` and :meth:`propagate` allows to
        estimate the "time discretization error". If the error is significant,
        a shorter time step shoud be used.

        Returns:
            qutip.solver.Result: Result of the propagation, using the same
            structure as :meth:`mesolve`.
        """
        if H is None:
            H = self.H
        if c_ops is None:
            c_ops = self.c_ops
        if e_ops is None:
            e_ops = []
        if args is None:
            args = {}
        result = QutipSolverResult()
        try:
            result.solver = propagator.__name__
        except AttributeError:
            try:
                result.solver = propagator.__class__.__name__
            except AttributeError:
                result.solver = 'n/a'
        result.times = np.array(tlist)
        result.states = []
        result.expect = []
        result.num_expect = len(e_ops)
        result.num_collapse = len(c_ops)
        for _ in e_ops:
            result.expect.append([])
        state = rho0
        if state is None:
            state = self.initial_state
        if len(e_ops) == 0:
            result.states.append(state)
        else:
            for (i, oper) in enumerate(e_ops):
                result.expect[i].append(expect(oper, state))
        controls = extract_controls([self])
        pulses_mapping = extract_controls_mapping([self], controls)
        mapping = pulses_mapping[0]  # "first objective" (dummy structure)
        pulses = [  # defined on the tlist intervals
            control_onto_interval(discretize(control, tlist, args=(args,)))
            for control in controls
        ]
        for time_index in range(len(tlist) - 1):  # index over intervals
            H_at_t = plug_in_pulse_values(H, pulses, mapping[0], time_index)
            c_ops_at_t = [
                plug_in_pulse_values(c_op, pulses, mapping[ic + 1], time_index)
                for (ic, c_op) in enumerate(c_ops)
            ]
            dt = tlist[time_index + 1] - tlist[time_index]
            state = propagator(
                H_at_t,
                state,
                dt,
                c_ops_at_t,
                initialize=True,  # initialize=(time_index == 0)
            )
            if len(e_ops) == 0:
                result.states.append(state)
            else:
                for (i, oper) in enumerate(e_ops):
                    result.expect[i].append(expect(oper, state))
        result.expect = [np.array(a) for a in result.expect]
        return result

    @classmethod
    def reset_symbol_counters(cls):
        """Reset the internal symbol counters used for printing objectives.

        See :meth:`summarize`.
        """
        cls._counter = defaultdict(int)
        cls._counter['u{count}(t)'] = 1
        cls._count_cache = {}

    def summarize(self, use_unicode=True, reset_symbol_counters=False):
        """Return a one-line summary of the objective as a string.

        Args:
            use_unicode(bool): If False, only use ascii symbols in the output
            reset_symbol_counters(bool): If True, reset the internal object
                counters (see :meth:`reset_symbol_counters`) before calculating
                the result

        The :meth:`summarize` method (which is also used for the :func:`repr`
        and ``__str__`` of an :class:`Objective`) keeps per-process internal
        counters for the various categories of objects that may occur as
        attributes of an :class:`Objective` (kets, bras, Hermitian operators,
        non-Hermitian Operators, density matrices, Liouvillians, Lindblad
        operators, numpy arrays, control functions). This allows to keep track
        of objects across multiple objectives. The counters can be reset with
        :meth:`reset_symbol_counters`.

        The ouput uses various unicode symbols (or ascii-equivalents, if
        `use_unicode` is False):

        - 'Ψ' ('Psi') for :class:`qutip.Qobj` quantum states (kets or bras)
        - 'ρ' ('rho') for :class:`qutip.Qobj` operators that occur as initial
          or target states (density matrices)
        - 'L' for Lindblad operators (elements of :attr:`c_ops`)
        - 'H' for Hermitian :class:`qutip.Qobj` operators (Hamiltonians)
        - 'A' for non-Hermitian :class:`qutip.Qobj` operators in :attr:`H`
        - '𝓛' ('Lv') for :class:`qutip.Qobj` super-operators (Liouvillians)
        - 'a' for numpy arrays (of any dimension)
        - 'u' for (callable) control functions.

        Example:

            >>> from qutip import ket, tensor, sigmaz, sigmax, sigmap, identity
            >>> u1 = lambda t, args: 1.0
            >>> u2 = lambda t, args: 1.0
            >>> a1 = np.random.random(100) + 1j*np.random.random(100)
            >>> a2 = np.random.random(100) + 1j*np.random.random(100)
            >>> H = [
            ...     tensor(sigmaz(), identity(2)) +
            ...     tensor(identity(2), sigmaz()),
            ...     [tensor(sigmax(), identity(2)), u1],
            ...     [tensor(identity(2), sigmax()), u2]]
            >>> C1 = [[tensor(identity(2), sigmap()), a1]]
            >>> C2 = [[tensor(sigmap(), identity(2)), a2]]
            >>> ket00 = ket((0,0))
            >>> ket11 = ket((1,1))
            >>> obj = Objective(
            ...     initial_state=ket00,
            ...     target=ket11,
            ...     H=H
            ... )
            >>> obj.reset_symbol_counters()
            >>> obj.summarize()
            '|Ψ₀(2⊗2)⟩ to |Ψ₁(2⊗2)⟩ via [H₀[2⊗2,2⊗2], [H₁[2⊗2,2⊗2], u₁(t)], [H₂[2⊗2,2⊗2], u₂(t)]]'
            >>> obj = Objective(
            ...     initial_state=ket00,
            ...     target=ket11,
            ...     H=H,
            ...     c_ops=[C1, C2]
            ... )
            >>> obj.summarize()
            '|Ψ₀(2⊗2)⟩ to |Ψ₁(2⊗2)⟩ via {H:[H₀[2⊗2,2⊗2], [H₁[2⊗2,2⊗2], u₁(t)], [H₂[2⊗2,2⊗2], u₂(t)]], c_ops:([[L₀[2⊗2,2⊗2], a₀[100]]],[[L₁[2⊗2,2⊗2], a₁[100]]])}'
            >>> obj.summarize(use_unicode=False)
            '|Psi0(2*2)> to |Psi1(2*2)> via {H:[H0[2*2,2*2], [H1[2*2,2*2], u1(t)], [H2[2*2,2*2], u2(t)]], c_ops:([[L0[2*2,2*2], a0[100]]],[[L1[2*2,2*2], a1[100]]])}'
            >>> copy.deepcopy(obj).summarize()  # different objects!
            '|Ψ₂(2⊗2)⟩ to |Ψ₃(2⊗2)⟩ via {H:[H₃[2⊗2,2⊗2], [H₄[2⊗2,2⊗2], u₁(t)], [H₅[2⊗2,2⊗2], u₂(t)]], c_ops:([[L₂[2⊗2,2⊗2], a₂[100]]],[[L₃[2⊗2,2⊗2], a₃[100]]])}'
            >>> copy.deepcopy(obj).summarize(reset_symbol_counters=True)
            '|Ψ₀(2⊗2)⟩ to |Ψ₁(2⊗2)⟩ via {H:[H₀[2⊗2,2⊗2], [H₁[2⊗2,2⊗2], u₁(t)], [H₂[2⊗2,2⊗2], u₂(t)]], c_ops:([[L₀[2⊗2,2⊗2], a₀[100]]],[[L₁[2⊗2,2⊗2], a₁[100]]])}'
        """
        if reset_symbol_counters:
            self.reset_symbol_counters()
        res = _summarize_component(
            self.initial_state,
            role='state',
            counter=self._counter,
            count_cache=self._count_cache,
            use_unicode=use_unicode,
        )
        if self.target is not None:
            res += " to "
            role = 'target'
            if (
                isinstance(self.initial_state, qutip.Qobj)
                and isinstance(self.target, qutip.Qobj)
                and (self.target.dims == self.initial_state.dims)
            ):
                role = 'state'
            res += _summarize_component(
                self.target,
                role=role,
                counter=self._counter,
                count_cache=self._count_cache,
                use_unicode=use_unicode,
            )
        res += " via "

        if len(self.c_ops) == 0:
            res += _summarize_component(
                self.H,
                role='op',
                counter=self._counter,
                count_cache=self._count_cache,
                use_unicode=use_unicode,
            )
        else:
            res += '{H:'
            res += _summarize_component(
                self.H,
                role='op',
                counter=self._counter,
                count_cache=self._count_cache,
                use_unicode=use_unicode,
            )
            rendered_c_ops = [
                _summarize_component(
                    c_op,
                    role='lindblad',
                    counter=self._counter,
                    count_cache=self._count_cache,
                    use_unicode=use_unicode,
                )
                for c_op in self.c_ops
            ]
            res += ', c_ops:('
            res += ",".join(rendered_c_ops)
            res += ')}'
        return res

    def __str__(self):
        return self.summarize(use_unicode=self.str_use_unicode)

    def __repr__(self):
        return "%s[%s]" % (self.__class__.__name__, str(self))


def _Objective_reduce_init(initial_state, H, target, c_ops):
    # args-only version of Objective.__init__, for _Objective_reduce
    return Objective(
        initial_state=initial_state, H=H, target=target, c_ops=c_ops
    )


def _Objective_reduce(obj):
    """Reduce :class:`Objective` for pickling.

    This is a reduction function for customized pickling, see
    :func:`copyreg.pickle`. It is used in :meth:`.Result.dump`.

    In the standard-library-pickle, lambdas are not pickleable, so we replace
    those non-pickleable entries with a placeholder.
    """
    return (
        _Objective_reduce_init,
        (
            obj.initial_state,
            _remove_functions_from_nested_list(obj.H),
            obj.target,
            _remove_functions_from_nested_list(obj.c_ops),
        ),
        {
            k: v
            for (k, v) in obj.__dict__.items()
            if k not in obj._default_attribs
        },
    )


class _ControlPlaceholder:
    """Placeholder for a control function, for pickling."""

    def __init__(self, id):
        self.id = id

    def __str__(self):
        return "u%s" % self.id

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.id)

    def __eq__(self, other):
        return (self.__class__ == other.__class__) and (self.id == other.id)


def _remove_functions_from_nested_list(lst):
    if isinstance(lst, list):
        return [_remove_functions_from_nested_list(v) for v in lst]
    else:
        if callable(lst) and not isinstance(lst, qutip.Qobj):
            return _ControlPlaceholder(id(lst))
        else:
            return lst


def _plug_in_array_controls_as_func(H, controls, mapping, tlist):
    """Convert array controls to piece-wise constant functions

    It uses the piece-wise constant convention of mesolve: pulses switch in the
    middle between to `tlist` points.

    This is a workaround for https://github.com/qutip/qutip/issues/932
    """
    H = _nested_list_shallow_copy(H)
    T = tlist[-1]
    nt = len(tlist)
    for (control, control_mapping) in zip(controls, mapping):
        if isinstance(control, np.ndarray):
            for i in control_mapping:
                # Use the same formula that QuTiP normally passes to Cython for
                # compilation
                H[i][1] = partial(_array_as_func, array=control, T=T, nt=nt)
        else:
            continue
    return H


def _array_as_func(t, args, array, T, nt):
    return (
        0
        if (t > float(T))
        else array[int(round(float(nt - 1) * (t / float(T))))]
    )


def _reversed_enumerate(collection):
    # better than `reversed(list(enumerate(collection)))`, because it doesn't
    # create copies, via `list`
    return zip(reversed(range(len(collection))), reversed(collection))


def _rho1(basis_states):
    """State ρ₁ from the "3states" functional"""
    d = len(basis_states)  # dimension of logical subspace
    return sum(
        [
            (2 * (d - i) / (d * (d + 1))) * psi * psi.dag()
            for (i, psi) in enumerate(basis_states)
            # note that i is 0-based, unlike in the paper
        ]
    )


def _rho2(basis_states):
    """State ρ₂ from the "3states" functional"""
    d = len(basis_states)  # dimension of logical subspace
    return (1.0 / d) * sum(
        [
            psi_i * psi_j.dag()
            for (psi_i, psi_j) in itertools.product(basis_states, repeat=2)
        ]
    )


def _rho3(basis_states):
    """State ρ₃ from the "3states" functional"""
    d = len(basis_states)  # dimension of logical subspace
    return (1.0 / d) * sum([psi * psi.dag() for psi in basis_states])


def gate_objectives(
    basis_states,
    gate,
    H,
    *,
    c_ops=None,
    local_invariants=False,
    liouville_states_set=None,
    weights=None,
    normalize_weights=True
):
    r"""Construct a list of objectives for optimizing towards a quantum gate

    Args:
        basis_states (list[qutip.Qobj]): A list of $n$ canonical basis states
        gate: The gate to optimize for, as a $n \times n$ matrix-like object
            (must have a `shape` attribute, and be indexable by two indices).
            Alternatively, `gate` may be the string 'perfect_entangler' or
            'PE', to indicate the optimization for an arbitrary two-qubit
            perfect entangler.
        H (list or qutip.Qobj): The Hamiltonian (or Liouvillian) for the time
            evolution, in nested-list format.
        c_ops (list or None): A list of collapse (Lindblad) operators, or None
            for unitary dynamics or if `H` is a Liouvillian (preferred!)
        local_invariants (bool): If True, initialize the objectives for an
            optimization towards a two-qubit gate that is "locally equivalent"
            to `gate`. That is, the result of the optimization should implement
            `gate` up to single-qubit operations.
        liouville_states_set (None or str): If not None, one of "full",
            "3states", "d+1". This sets the objectives for a gate
            optimization in Liouville space, using the states defined in
            Goerz et al. New J. Phys. 16, 055012 (2014). See Examples for
            details.
        weights (None or list): If given as a list, weights for the different
            objectives. These will be added as a custom attribute to the
            respective :class:`.Objective`, and may be used by a particular
            functional (`chi_constructor`). The intended use case is for the
            `liouville_states_set` values '3states', and 'd+1', where the
            different objectives have clear physical interpretations that might
            be given differing importance. A weight of 0 will completely drop
            the corresponding objective.
        normalize_weights (bool): If True, and if `weights` is given as a list
            of values, normalize the weights so that they sum to $N$, the
            number of objectives. IF False, the weights will be used unchanged.

    Returns:
        list[Objective]: The objectives that define the optimization towards
        the gate. For a "normal" gate with a basis in Hilbert space, the
        objectives will have the `basis_states` as each
        :attr:`~.Objective.initial_state` and the result of applying `gate`
        to the `basis_states` as each :attr:`~.Objective.target`.

        For an optimization towards a perfect-entangler, or for the
        `local_invariants` of the given `gate`, each
        :attr:`~.Objective.initial_state` will be the Bell basis state
        described in "Theorem 1" in Y. Makhlin, Quantum Inf. Process. 1, 243
        (2002), derived from the canonical `basis_states`. The
        :attr:`~.Objective.target` will be the string 'PE' for a
        perfect-entanglers optimization, and `gate` for the local-invariants
        optimization.

    Raises:
        ValueError: If `gate`, `basis_states`, and `local_invariants` are
            incompatible, or `gate` is invalid (not a recognized string)

    .. Note::

        The dimension of the `basis_states` is not required to be the dimension
        of the `gate`; the `basis_states` may define a logical subspace in a
        larger Hilbert space.

    Examples:

        * A single-qubit gate::

            >>> from qutip import ket, bra, tensor
            >>> from qutip import sigmaz, sigmax, sigmay, sigmam, identity
            >>> basis = [ket([0]), ket([1])]
            >>> gate = sigmay()  # = -i|0⟩⟨1| + i|1⟩⟨0|
            >>> H = [sigmaz(),[sigmax(), lambda t, args: 1.0]]
            >>> objectives = gate_objectives(basis, gate, H)
            >>> assert objectives == [
            ...     Objective(
            ...         initial_state=basis[0],
            ...         target=(1j * basis[1]),
            ...         H=H
            ...     ),
            ...     Objective(
            ...         initial_state=basis[1],
            ...         target=(-1j * basis[0]),
            ...         H=H
            ...     )
            ... ]

        * An arbitrary two-qubit perfect entangler:

            >>> basis = [ket(n) for n in [(0, 0), (0, 1), (1, 0), (1, 1)]]
            >>> H = [
            ...     tensor(sigmaz(), identity(2)) +
            ...     tensor(identity(2), sigmaz()),
            ...     [tensor(sigmax(), identity(2)), lambda t, args: 1.0],
            ...     [tensor(identity(2), sigmax()), lambda t, args: 1.0]]
            >>> objectives = gate_objectives(basis, 'PE', H)
            >>> from weylchamber import bell_basis
            >>> for i in range(4):
            ...     assert objectives[i] == Objective(
            ...        initial_state=bell_basis(basis)[i],
            ...        target='PE',
            ...        H=H
            ...     )

        * A two-qubit gate, up to single-qubit operation ("local invariants"):

            >>> CNOT = qutip.Qobj(
            ...     [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
            ...     dims=[[2, 2], [2, 2]]
            ... )
            >>> objectives = gate_objectives(
            ...     basis, CNOT, H, local_invariants=True
            ... )
            >>> for i in range(4):
            ...     assert objectives[i] == Objective(
            ...        initial_state=bell_basis(basis)[i],
            ...        target=CNOT,
            ...        H=H
            ...     )

        * A two-qubit gate in a dissipative system tracked by 3 density
          matrices::

            >>> L = krotov.objectives.liouvillian(H, c_ops=[
            ...     tensor(sigmam(), identity(2)),
            ...     tensor(identity(2), sigmam())])
            >>> objectives = gate_objectives(
            ...     basis, CNOT, L,
            ...     liouville_states_set='3states',
            ...     weights=[20, 1, 1]
            ... )

          The three states, for a system with a logical subspace of dimension
          $d$ with a basis $\{\ket{i}\}$, $i \in [1, d]$ are:

          .. math::

            \Op{\rho}_1 &= \sum_{i=1}^{d}
                \frac{2 (d-i+1)}{d (d+1)} \ketbra{i}{i} \\
            \Op{\rho}_2 &= \sum_{i,j=1}^{d}
                \frac{1}{d} \ketbra{i}{j} \\
            \Op{\rho}_3 &= \sum_{i=1}^{d}
                \frac{1}{d} \ketbra{i}{i}

          The explicit form of the three states in this example is::

            >>> assert np.allclose(objectives[0].initial_state.full(),
            ...     np.diag([0.4, 0.3, 0.2, 0.1]))

            >>> assert np.allclose(objectives[1].initial_state.full(),
            ...     np.full((4, 4), 1/4))

            >>> assert np.allclose(objectives[2].initial_state.full(),
            ...     np.diag([1/4, 1/4, 1/4, 1/4]))

          The objectives in this example are weighted (20/1/1)::

            >>> "%.5f" % objectives[0].weight
            '2.72727'
            >>> "%.5f" % objectives[1].weight
            '0.13636'
            >>> "%.5f" % objectives[2].weight
            '0.13636'
            >>> sum_of_weights = sum([obj.weight for obj in objectives])
            >>> "%.1f" % sum_of_weights
            '3.0'

        * A two-qubit gate in a dissipative system tracked by $d + 1 = 5$
          pure-state density matrices::

            >>> objectives = gate_objectives(
            ...     basis, CNOT, L,
            ...     liouville_states_set='d+1'
            ... )

          The first four `initial_states` are the pure states corresponding to
          the Hilbert space basis

            >>> assert objectives[0].initial_state == qutip.ket2dm(ket('00'))
            >>> assert objectives[1].initial_state == qutip.ket2dm(ket('01'))
            >>> assert objectives[2].initial_state == qutip.ket2dm(ket('10'))
            >>> assert objectives[3].initial_state == qutip.ket2dm(ket('11'))

          The fifth state is $\Op{\rho}_2$ from '3states'::

            >>> assert np.allclose(objectives[4].initial_state.full(),
            ...     np.full((4, 4), 1/4))

        * A two-qubit gate in a dissipative system tracked by the full
          Liouville space basis::

            >>> objectives = gate_objectives(
            ...     basis, CNOT, L,
            ...     liouville_states_set='full'
            ... )

          The Liouville space basis states are all the possible dyadic products
          of the Hilbert space basis::

            >>> assert objectives[0].initial_state == ket('00') * bra('00')
            >>> assert objectives[1].initial_state == ket('00') * bra('01')
            >>> assert objectives[2].initial_state == ket('00') * bra('10')
            >>> assert objectives[3].initial_state == ket('00') * bra('11')
            >>> assert objectives[4].initial_state == ket('01') * bra('00')
            >>> assert objectives[5].initial_state == ket('01') * bra('01')
            >>> assert objectives[6].initial_state == ket('01') * bra('10')
            >>> assert objectives[7].initial_state == ket('01') * bra('11')
            >>> assert objectives[8].initial_state == ket('10') * bra('00')
            >>> assert objectives[9].initial_state == ket('10') * bra('01')
            >>> assert objectives[10].initial_state == ket('10') * bra('10')
            >>> assert objectives[11].initial_state == ket('10') * bra('11')
            >>> assert objectives[12].initial_state == ket('11') * bra('00')
            >>> assert objectives[13].initial_state == ket('11') * bra('01')
            >>> assert objectives[14].initial_state == ket('11') * bra('10')
            >>> assert objectives[15].initial_state == ket('11') * bra('11')
    """
    if isinstance(gate, str):
        if gate.lower().replace(' ', '_') in ['pe', 'perfect_entangler']:
            return _gate_objectives_li_pe(basis_states, 'PE', H, c_ops)
        else:
            raise ValueError(
                "gate must be either a square matrix, or one of the strings "
                "'PE' or 'perfect_entangler', not '" + gate + "'"
            )
    elif local_invariants:
        if not gate.shape == (4, 4):
            raise ValueError(
                "If local_invariants is True, gate must be a 4 × 4 matrix, "
                "not " + str(gate.shape)
            )
        return _gate_objectives_li_pe(basis_states, gate, H, c_ops)

    # "Normal" gate:

    if not gate.shape[0] == gate.shape[1] == len(basis_states):
        raise ValueError(
            "gate must be a matrix of the same dimension as the number of "
            "basis states"
        )
    mapped_basis_states = [
        sum(
            [
                complex(gate[i, j]) * basis_states[i]
                for i in range(gate.shape[0])
            ]
        )
        for j in range(gate.shape[1])
    ]
    # Lots of gates just rearrange the basis states, and we can avoid some
    # complexity (and make the repr of an Objective look nicer) by identifying
    # this and setting the mapped_basis_states to the identical objects as the
    # original basis_states
    for (i, state) in enumerate(mapped_basis_states):
        for (j, basis_state) in enumerate(basis_states):
            if state == basis_state:
                mapped_basis_states[i] = basis_state
    if liouville_states_set is None:
        # standard gate in Hilbert space
        initial_states = basis_states
        target_states = mapped_basis_states
    elif liouville_states_set.lower() == 'full':
        initial_states = [
            psi_i * psi_j.dag()
            for (psi_i, psi_j) in itertools.product(basis_states, repeat=2)
        ]
        target_states = [
            psi_i * psi_j.dag()
            for (psi_i, psi_j) in itertools.product(
                mapped_basis_states, repeat=2
            )
        ]
    elif liouville_states_set.replace(" ", "").lower() == '3states':
        d = len(basis_states)  # dimension of logical subspace
        initial_states = [
            _rho1(basis_states),
            _rho2(basis_states),
            _rho3(basis_states),
        ]
        target_states = [
            _rho1(mapped_basis_states),
            _rho2(mapped_basis_states),
            _rho3(mapped_basis_states),
        ]
    elif liouville_states_set.replace(" ", "").lower() == 'd+1':
        d = len(basis_states)
        initial_states = [
            basis_states[i] * basis_states[i].dag() for i in range(d)
        ]
        initial_states.append(_rho2(basis_states))
        target_states = [
            mapped_basis_states[i] * mapped_basis_states[i].dag()
            for i in range(d)
        ]
        target_states.append(_rho2(mapped_basis_states))
    else:
        raise ValueError(
            "Invalid `liouville_states_set`: %s" % liouville_states_set
        )
    objectives = [
        Objective(
            initial_state=initial_state, target=target_state, H=H, c_ops=c_ops
        )
        for (initial_state, target_state) in zip(initial_states, target_states)
    ]
    # apply weights
    if weights is not None:
        if len(weights) != len(objectives):
            raise ValueError(
                "If weight are given, there must be a weight for each "
                "objective"
            )
        if normalize_weights:
            N = len(objectives)
            weights = N * np.array(weights) / np.sum(weights)
        for (i, weight) in _reversed_enumerate(weights):
            weight = float(weight)
            if weight < 0:
                raise ValueError("weights must be greater than zero")
            objectives[i].weight = weight
            if weight == 0:
                del objectives[i]
    return objectives


def _gate_objectives_li_pe(basis_states, gate, H, c_ops):
    """Objectives for two-qubit local-invariants or perfect-entangler
    optimizaton"""
    if len(basis_states) != 4:
        raise ValueError(
            "Optimization towards a two-qubit gate requires 4 basis_states"
        )
    # Bell states as in "Theorem 1" in
    # Y. Makhlin, Quantum Inf. Process. 1, 243 (2002)
    psi1 = (basis_states[0] + basis_states[3]) / np.sqrt(2)
    psi2 = (1j * basis_states[1] + 1j * basis_states[2]) / np.sqrt(2)
    psi3 = (basis_states[1] - basis_states[2]) / np.sqrt(2)
    psi4 = (1j * basis_states[0] - 1j * basis_states[3]) / np.sqrt(2)
    return [
        Objective(initial_state=psi, target=gate, H=H, c_ops=c_ops)
        for psi in [psi1, psi2, psi3, psi4]
    ]


def ensemble_objectives(objectives, Hs, *, keep_original_objectives=True):
    """Extend `objectives` for an "ensemble optimization"

    This creates a list of objectives for an optimization for robustness with
    respect to variations in some parameter of the Hamiltonian. The trick is to
    simply optimize over the average of multiple copies of the system
    (the `Hs`) sampling that variation. See
    Goerz, Halperin, Aytac, Koch, Whaley. Phys. Rev. A 90, 032329 (2014)
    for details.

    Args:
        objectives (list[Objective]): The $n$ original objectives
        Hs (list): List of $m$ variations of the original
            Hamiltonian/Liouvillian
        keep_original_objectives (bool): If given as False, drop the original
            objectives from the result. This is especially useful if `Hs`
            contains the original Hamiltonian (which is often more
            straightforward)

    Returns:
        list[Objective]: List of $n (m+1)$ new objectives that consists of the
        original objectives, plus one copy of the original objectives per
        element of `Hs` where the `H` attribute of each objectives is
        replaced by that element. Alternatively, for
        ``keep_original_objectives=False``, list of $n m$ new objectives
        without the original objectives.
    """
    new_objectives = []
    if keep_original_objectives:
        new_objectives = copy.copy(objectives)
    for H in Hs:
        for obj in objectives:
            new_objectives.append(
                Objective(
                    H=H,
                    initial_state=obj.initial_state,
                    target=obj.target,
                    c_ops=obj.c_ops,
                )
            )
    return new_objectives


def liouvillian(H, c_ops):
    """Convert Hamiltonian and Lindblad operators into a Liouvillian.

    This is like :func:`qutip.superoperator.liouvillian`, but `H` may be a
    time-dependent Hamiltonian in nested-list format. `H` is assumed to contain
    a drift Hamiltonian, and the Lindblad operators in `c_ops` cannot be
    time-dependent.
    """
    if isinstance(H, qutip.Qobj):
        return qutip.liouvillian(H, c_ops)
    elif isinstance(H, list):
        res = []
        for spec in H:
            if isinstance(spec, qutip.Qobj):
                res.append(qutip.liouvillian(spec, c_ops))
                c_ops = []
            else:
                res.append([qutip.liouvillian(spec[0]), spec[1]])
        assert len(c_ops) == 0, "No drift Hamiltonian"
        return res
    else:
        raise ValueError(
            "H must either be a Qobj, or a time-dependent Hamiltonian in "
            "nested-list format"
        )


def _obj_str_pattern(obj, role, use_unicode):
    """For a given `obj` return a "string pattern" to be used summerizing that
    `obj` as the component of an Objective.

    The resulting string should be formatted with the keys ``count`` (a unique
    str counter for the `obj` in its category), and ``dims`` with a str summary
    of the shape/dimensions of `obj`.
    """
    if callable(obj) and not isinstance(obj, qutip.Qobj):
        if role == 'op':
            return 'u{count}(t)'  # control
        else:
            return None
    elif isinstance(obj, np.ndarray):
        return 'a{count}[{dims}]'  # array
    elif isinstance(obj, qutip.Qobj):
        if obj.type == 'ket':
            if use_unicode:
                return '|Ψ{count}({dims})⟩'
            else:
                return '|Psi{count}({dims})>'
        elif obj.type == 'bra':
            if use_unicode:
                return '⟨Ψ{count}({dims})|'
            else:
                return '<Psi{count}({dims})|'
        elif obj.type == 'oper':
            if role == 'lindblad':
                return 'L{count}[{dims}]'
            if obj.isherm:
                if role == 'state':
                    if use_unicode:
                        return 'ρ{count}[{dims}]'
                    else:
                        return 'rho{count}[{dims}]'
                else:
                    return 'H{count}[{dims}]'
            else:
                return 'A{count}[{dims}]'
        elif obj.type == 'super':
            if use_unicode:
                return "𝓛{count}[{dims}]"
            else:
                return "Lv{count}[{dims}]"
        else:
            raise NotImplementedError("Unknown qobj type: %s" % obj.type)
    else:
        return None  # unknown object


def _obj_dims_str(obj, use_unicode):
    """Return a string that summarizes the shape/dimensions of the given obj"""
    tensor = '*'
    if use_unicode:
        tensor = '⊗'
    if isinstance(obj, qutip.Qobj):
        if obj.type == 'ket':
            return tensor.join(["%d" % d for d in obj.dims[0]])
        elif obj.type == 'bra':
            return tensor.join(["%d" % d for d in obj.dims[1]])
        elif obj.type == 'oper':
            return ",".join(
                [tensor.join(["%d" % d for d in dim]) for dim in obj.dims]
            )
        elif obj.type == 'super':
            dims = []
            for dim in obj.dims:
                dim1 = tensor.join(["%d" % d for d in dim[0]])
                dim2 = tensor.join(["%d" % d for d in dim[1]])
                dims.append('[%s,%s]' % (dim1, dim2))
            return ",".join(dims)
        else:
            raise NotImplementedError("Unknown qobj type: %s" % obj.type)
    elif hasattr(obj, 'shape'):
        return ",".join(str(int(d)) for d in obj.shape)
    else:
        return None


def _summarize_component(
    obj, role, counter=None, count_cache=None, use_unicode=True
):
    """Return a string for an `obj` that appears as a component when
    summarizing an Objective.


    Args:
        obj: Any object that occurs in any attribute of an :class:`Objective`.
        role (str): The role the `obj` plays in the objective, one of 'state',
            'target', 'op', 'lindblad'.
        counter (defaultdict): A dictionary that matches a string pattern (as
            an identifier for the category of the `obj` to the number of
            objects that have been observed in that category.
        count_cache (dict): For any previously seen `obj`, the value of
            `counter[_obj_str_pattern(obj)]` from when that `obj` was first
            seen
        use_unicode (bool): If true, return a unicode representation, ascii
            otherwise

    For a given `obj`, let's say a ket-:class:`qutip.Qobj` this works as
    follows:

    - :func:`_obj_str_pattern` returns ``"|Ψ{count}({dims})⟩"``. This
      unformatted string defines the "category" of the `obj`. Any object with
      the same unformatted string is in the same category (kets, bras,
      Hermitian operators, Non-Hermitian Operators, density matrices,
      Liouvillian, Lindblad operators)
    - if we have seen that particular `obj` before, ``{count}`` gets set from
      the `count_cache`. Otherwise, we increase the count for the category, and
      add the object to `count_cache` (using the memory-address of `obj` as a
      key)
    - the ``{dims}}`` get set by :func:`_obj_dims_str`.
    """
    allowed_roles = ['state', 'target', 'op', 'lindblad']
    if role not in allowed_roles:
        raise ValueError("Unknown %s not in %s" % (role, allowed_roles))
    if counter is None:
        counter = Objective._counter
    if count_cache is None:
        count_cache = Objective._count_cache
    if isinstance(obj, list):
        return _summarize_nested_list(
            obj, role, counter, count_cache, use_unicode
        )
    elif isinstance(obj, _ControlPlaceholder):
        return str(obj)
    elif isinstance(obj, (float, complex)):
        return str(obj)
    str_pattern = _obj_str_pattern(obj, role, use_unicode)
    if str_pattern is None:
        # we're dealing with some kind of unknown object, so we'll fall back to
        # a single line, truncated str-representation
        res = str(obj).replace("\n", " ")
        if len(res) > 40:
            if use_unicode:
                res = res[:39] + "…"
            else:
                res = res[:37] + "..."
    else:
        key = id(obj)
        # Lot's of objectes that occur as components of an Objectives, e.g.
        # numpy arrays, are not hashable, and thus cannot be used as dictionary
        # keys. Thus we use the memory address of the object (integer returned
        # by :func:`id` as the key. This matches the intuition that objects
        # with the same subscript are really the same object, not different
        # objects that are potentially numerically equal.
        dims = _obj_dims_str(obj, use_unicode)
        if key in count_cache:
            count = count_cache[key]
        else:
            count = counter[str_pattern]
            count_cache[key] = count
            counter[str_pattern] += 1
            # Things look a bit nicer if we synchronize the counters for
            # Hermitian and Non-Hermitian operators
            if str_pattern == 'A{count}[{dims}]':
                counter['H{count}[{dims}]'] += 1
            elif str_pattern == 'H{count}[{dims}]':
                counter['A{count}[{dims}]'] += 1
        count_str = str(count)
        if use_unicode:
            # transform all digits in counter to unicode subscripts.
            # Subscript symbols start at code point 0x2080
            count_str = "".join(
                chr(ord(d) - ord('0') + 0x2080) for d in str(count)
            )
        res = str_pattern.format(count=count_str, dims=dims)
    return res


def _summarize_nested_list(lst, role, counter, count_cache, use_unicode):
    """Recursively call :func:`summarize_components` for the nested lists that
    qutip uses for time-dependent operators
    """
    return (
        '['
        + ", ".join(
            [
                _summarize_component(
                    obj, role, counter, count_cache, use_unicode
                )
                for obj in lst
            ]
        )
        + ']'
    )


def _recursive_eq(a, b):
    """Recursively compare `a` and `b`.

    The parameters `a` and `b` are assumed to be attributes of
    :class:`Objective`, so we're making some assumptions about these either
    being "standard" values (:class:`qutip.Qobj` and numpy arrays potentially
    occuring in nested lists as time-dependent controls, or else `a` and `b`
    comparing directly with ``==`` without throwing an exception.

    If that's not enough, it will be up to the user to use wrapper
    objects that implement a useful ``__eq__``.

    We need this because '==' does not work for numpy arrays.
    Cf. https://stackoverflow.com/questions/55778847
    """
    if type(a) != type(b):
        return False
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(
            _recursive_eq(v, w) for (v, w) in zip(a, b)
        )
    elif isinstance(a, dict):
        return len(a) == len(b) and all(
            _recursive_eq(v, b[k]) for (k, v) in a.items()
        )
    else:
        try:
            return bool(a == b)
        except ValueError:
            if isinstance(a, np.ndarray):
                return np.array_equal(a, b)
            else:
                raise
