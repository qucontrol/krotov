import sys
import copy
from functools import partial

import qutip
from qutip.solver import Result as QutipSolverResult
import numpy as np

from .structural_conversions import (
    _nested_list_shallow_copy, extract_controls, extract_controls_mapping,
    plug_in_pulse_values, control_onto_interval, discretize)

__all__ = [
    'Objective', 'summarize_qobj', 'CtrlCounter', 'gate_objectives',
    'ensemble_objectives']


FIX_QUTIP_932 = True
"""Workaround for `QuTiP issue 932`_.

If True, and only when running on macOS, in :meth:`Objective.mesolve`,
replace any array controls with an equivalent function. This results in a
signficant slowdown of the propagation, as it circumvents the use of Cython.

.. _QuTiP issue 932: https://github.com/qutip/qutip/issues/932
"""


def _adjoint(op):
    """Calculate adjoint of an operator in QuTiP nested-list format.
    Controls are not modified."""
    if isinstance(op, list):
        adjoint_op = []
        for item in op:
            if isinstance(item, list):
                assert len(item) == 2
                adjoint_op.append([item[0].dag(), item[1]])
            else:
                adjoint_op.append(item.dag())
        return adjoint_op
    elif op is None:
        return None
    else:
        return op.dag()


def CtrlCounter():
    """Constructor for a counter of controls.

    Returns a callable that returns a unique integer (starting at 1) for every
    distinct control that is passed to it. This is intended for use with
    :func:`summarize_qobj`.

    Example:

        >>> ctrl_counter = CtrlCounter()
        >>> ctrl1 = np.zeros(10)
        >>> ctrl_counter(ctrl1)
        1
        >>> ctrl2 = np.zeros(10)
        >>> assert ctrl2 is not ctrl1
        >>> ctrl_counter(ctrl2)
        2
        >>> ctrl_counter(ctrl1)
        1
        >>> ctrl3 = lambda t, args: 0.0
        >>> ctrl_counter(ctrl3)
        3
    """

    counter = 0
    registry = {}

    def ctrl_counter(ctrl):
        nonlocal counter
        if isinstance(ctrl, np.ndarray):
            ctrl = id(ctrl)
        if ctrl not in registry:
            counter += 1
            registry[ctrl] = counter
        return registry[ctrl]

    return ctrl_counter


_CTRL_COUNTER = CtrlCounter()  #: internal counter for controls


class Objective:
    """A single objective for optimization with Krotov's method

    Args:
        initial_state (qutip.Qobj): value for :attr:`initial_state`
        H (qutip.Qobj or list): value for :attr:`H`
        target (qutip.Qobj or None): value for :attr:`target`
        c_ops (list or None): value for :attr:`c_ops`

    Attributes:
        H (qutip.Qobj or list): The (time-dependent) Hamiltonian,
            cf. :func:`qutip.mesolve.mesolve`. This includes the control
            fields.
        initial_state (qutip.Qobj): The initial state
        target: An object describing the "target" of the optimization, for the
            dynamics starting from :attr:`initial_state`. Usually, this will be
            the target state (the state into which :attr:`initial_state` should
            evolve). More generally, it can be an arbitrary object meeting the
            requirements of a specific `chi_constructor` function that will be
            passed to :func:`.optimize_pulses`.
        c_ops (list or None): List of collapse operators,
            cf. :func:`~qutip.mesolve.mesolve`.

    Raises:
        ValueError: If any arguments have an invalid type
    """

    def __init__(self, *, initial_state, H, target, c_ops=None):
        if c_ops is None:
            c_ops = []
        if not isinstance(H, (qutip.Qobj, list)):
            raise ValueError(
                "Invalid H, must be a Qobj, or a nested list, not %s"
                % H.__class__.__name__
            )
        self.H = H
        if not isinstance(initial_state, qutip.Qobj):
            raise ValueError(
                "Invalid initial_state: must be Qobj, not %s"
                % initial_state.__class__.__name__
            )
        self.initial_state = initial_state
        self.target = target
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
        return Objective(
            H=_nested_list_shallow_copy(self.H),
            initial_state=self.initial_state,
            target=self.target,
            c_ops=[_nested_list_shallow_copy(c) for c in self.c_ops],
        )

    def __eq__(self, other):
        if other.__class__ is self.__class__:
            return (self.H, self.initial_state, self.target, self.c_ops) == (
                other.H,
                other.initial_state,
                other.target,
                other.c_ops,
            )
        else:
            return NotImplemented

    def __ne__(self, other):  # pragma: nocover
        result = self.__eq__(other)
        if result is NotImplemented:
            return NotImplemented
        else:
            return not result

    @property
    def adjoint(self):
        """The :class:`Objective` containing the adjoint of all components.

        This does not affect the controls in :attr:`H`: these are
        assumed to be real-valued. Also, :attr:`.Objective.target` will be left
        unchanged if it is not a :class:`qutip.Qobj`.
        """
        return Objective(
            H=_adjoint(self.H),
            initial_state=_adjoint(self.initial_state),
            target=(
                _adjoint(self.target)
                if isinstance(self.target, qutip.Qobj)
                else self.target
            ),
            c_ops=[_adjoint(op) for op in self.c_ops],
        )

    def mesolve(self, tlist, rho0=None, e_ops=None, **kwargs):
        """Run :func:`qutip.mesolve.mesolve` on the system of the objective

        Solve the dynamics for the :attr:`H` and :attr:`c_ops` of the
        objective. If `rho0` is not given, the :attr:`initial_state` will be
        propagated. All other arguments will be passed to
        :func:`qutip.mesolve.mesolve`.

        Returns:
            qutip.solver.Result: Result of the propagation, see
            :func:`qutip.mesolve.mesolve` for details.
        """
        if rho0 is None:
            rho0 = self.initial_state
        if e_ops is None:
            e_ops = []
        H = self.H
        c_ops = self.c_ops
        if FIX_QUTIP_932 and sys.platform == "darwin":  # pragma: no cover
            # "darwin" = macOS; the "pragma" excludes from coverage analysis
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
        return qutip.mesolve(
            H=H, rho0=rho0, tlist=tlist, c_ops=c_ops, e_ops=e_ops, **kwargs
        )

    def propagate(self, tlist, *, propagator, rho0=None, e_ops=None):
        """Propagates the system of the objective over the entire time grid

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
        in `tlist`)

        Comparing the result of :meth:`mesolve` and :meth:`propagate` allows to
        estimate the "time discretization error". If the error is significant,
        a shorter time step shoud be used.

        Returns:
            qutip.solver.Result: Result of the propagation, using the same
            structure as :meth:`mesolve`.
        """
        if e_ops is None:
            e_ops = []
        result = QutipSolverResult()
        result.solver = propagator.__name__
        result.times = copy.copy(tlist)
        result.states = []
        result.expect = []
        result.num_expect = len(e_ops)
        result.num_collapse = len(self.c_ops)
        for _ in e_ops:
            result.expect.append([])
        state = rho0
        if state is None:
            state = self.initial_state
        if len(e_ops) == 0:
            result.states.append(state)
        else:
            for (i, oper) in enumerate(e_ops):
                result.expect[i].append(qutip.expect(oper, state))
        controls = extract_controls([self])
        pulses_mapping = extract_controls_mapping([self], controls)
        mapping = pulses_mapping[0]  # "first objective" (dummy structure)
        pulses = [  # defined on the tlist intervals
            control_onto_interval(discretize(control, tlist))
            for control in controls
        ]
        for time_index in range(len(tlist) - 1):  # index over intervals
            H = plug_in_pulse_values(self.H, pulses, mapping[0], time_index)
            c_ops = [
                plug_in_pulse_values(c_op, pulses, mapping[ic + 1], time_index)
                for (ic, c_op) in enumerate(self.c_ops)
            ]
            dt = tlist[time_index + 1] - tlist[time_index]
            state = propagator(H, state, dt, c_ops)
            if len(e_ops) == 0:
                result.states.append(state)
            else:
                for (i, oper) in enumerate(e_ops):
                    result.expect[i].append(qutip.expect(oper, state))
        return result

    def summarize(self, ctrl_counter=None):
        """Return a one-line summary of the objective as a string

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
            >>> C1 = [tensor(identity(2), sigmap()), a1]
            >>> C2 = [tensor(sigmap(), identity(2)), a2]
            >>> ket00 = ket((0,0))
            >>> ket11 = ket((1,1))
            >>> obj = Objective(
            ...     initial_state=ket00,
            ...     target=ket11,
            ...     H=H
            ... )
            >>> obj.summarize()
            '|(2⊗2)⟩ - {[Herm[2⊗2,2⊗2], [Herm[2⊗2,2⊗2], u1(t)], [Herm[2⊗2,2⊗2], u2(t)]]} - |(2⊗2)⟩'
            >>> obj = Objective(
            ...     initial_state=ket00,
            ...     target=ket11,
            ...     H=H,
            ...     c_ops=[C1, C2]
            ... )
            >>> obj.summarize()
            '|(2⊗2)⟩ - {H:[Herm[2⊗2,2⊗2], [Herm[2⊗2,2⊗2], u1(t)], [Herm[2⊗2,2⊗2], u2(t)]], c_ops:([NonHerm[2⊗2,2⊗2], u3[complex128]],[NonHerm[2⊗2,2⊗2], u4[complex128]])} - |(2⊗2)⟩'
        """
        if ctrl_counter is None:
            ctrl_counter = _CTRL_COUNTER
        if len(self.c_ops) == 0:
            res = "%s - {%s}" % (
                summarize_qobj(self.initial_state, ctrl_counter),
                summarize_qobj(self.H, ctrl_counter),
            )
        else:
            res = "%s - {H:%s, c_ops:(%s)}" % (
                summarize_qobj(self.initial_state, ctrl_counter),
                summarize_qobj(self.H, ctrl_counter),
                ",".join(
                    [summarize_qobj(c_op, ctrl_counter) for c_op in self.c_ops]
                ),
            )
        if self.target is not None:
            if isinstance(self.target, qutip.Qobj):
                res += " - %s" % (summarize_qobj(self.target, ctrl_counter))
            else:
                res += " - %r" % self.target
        return res

    def __str__(self):
        return self.summarize()

    def __repr__(self):
        return "%s[%s]" % (self.__class__.__name__, self.summarize())

    def __getstate__(self):
        """Return data for the pickle serialization of an objective.

        This may not preserve time-dependent controls, and is only to enable
        the serialization of :class:`.Result` objects.
        """
        state = copy.copy(self.__dict__)
        # Remove the unpicklable entries.
        state['H'] = _remove_functions_from_nested_list(state['H'])
        state['c_ops'] = _remove_functions_from_nested_list(state['c_ops'])
        return state


class _ControlPlaceholder():
    """Placeholder for a control function, for pickling"""

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
        0 if (t > float(T)) else
        array[int(round(float(nt-1) * (t/float(T))))])


def summarize_qobj(obj, ctrl_counter=None):
    """Summarize a quantum object

    A counter created by :func:`CtrlCounter` may be passed to distinguish
    control fields.

    Example:

        >>> ket = qutip.ket([1, 0, 1])
        >>> summarize_qobj(ket)
        '|(2⊗2⊗2)⟩'
        >>> bra = ket.dag()
        >>> summarize_qobj(bra)
        '⟨(2⊗2⊗2)|'
        >>> rho = ket * bra
        >>> summarize_qobj(rho)
        'Herm[2⊗2⊗2,2⊗2⊗2]'
        >>> a = qutip.create(10)
        >>> summarize_qobj(a)
        'NonHerm[10,10]'
        >>> S = qutip.to_super(a)
        >>> summarize_qobj(S)
        '[[10,10],[10,10]]'
    """
    if ctrl_counter is None:
        ctrl_counter = _CTRL_COUNTER
    if isinstance(obj, list):
        return _summarize_qobj_nested_list(obj, ctrl_counter)
    elif callable(obj) and not isinstance(obj, qutip.Qobj):
        return 'u%d(t)' % ctrl_counter(obj)
    elif isinstance(obj, np.ndarray):
        return 'u%d[%s]' % (ctrl_counter(obj), obj.dtype.name)
    elif isinstance(obj, _ControlPlaceholder):
        return str(obj)
    elif not isinstance(obj, qutip.Qobj):
        raise TypeError("obj must be a Qobj, not %s" % obj.__class__.__name__)
    if obj.type == 'ket':
        dims = "⊗".join(["%d" % d for d in obj.dims[0]])
        return '|(%s)⟩' % dims
    elif obj.type == 'bra':
        dims = "⊗".join(["%d" % d for d in obj.dims[1]])
        return '⟨(%s)|' % dims
    elif obj.type == 'oper':
        dim1 = "⊗".join(["%d" % d for d in obj.dims[0]])
        dim2 = "⊗".join(["%d" % d for d in obj.dims[1]])
        if obj.isherm:
            return 'Herm[%s,%s]' % (dim1, dim2)
        else:
            return 'NonHerm[%s,%s]' % (dim1, dim2)
    elif obj.type == 'super':
        dims = []
        for dim in obj.dims:
            dim1 = "⊗".join(["%d" % d for d in dim[0]])
            dim2 = "⊗".join(["%d" % d for d in dim[1]])
            dims.append('[%s,%s]' % (dim1, dim2))
        return '[' + ",".join(dims) + ']'
    else:
        raise NotImplementedError("Unknown qobj type: %s" % obj.type)


def _summarize_qobj_nested_list(lst, ctrl_counter):
    """Summarize a nested-list time-dependent quantum object"""
    return (
        '[' +
        ", ".join([summarize_qobj(obj, ctrl_counter) for obj in lst]) +
        ']')


def gate_objectives(basis_states, gate, H, c_ops=None, local_invariants=False):
    r"""Construct a list of objectives for optimizing towards a quantum gate

    Args:
        basis_states (list[qutip.Qobj]): A list of $n$ canonical basis states
        gate: The gate to optimize for, as a $n \times n$ matrix-like object
            (must have a `shape` attribute, and be indexable by two indices).
            Alternatively, `gate` may be the string 'perfect_entangler' or
            'PE', to indicate the optimization for an arbitrary two-qubit
            perfect entangler.
        H (list or qutip.Qobj): The Hamiltonian for the time evolution
        c_ops (list or None): A list of collapse (Lindblad) operators, or None
            for unitary dynamics
        local_invariants (bool): If True, initialize the objectives for an
            optimization towards a two-qubit gate that is "locally equivalent"
            to `gate`. That is, the result of the optimization should implement
            `gate` up to single-qubit operations.

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

            >>> from qutip import ket, tensor, sigmaz, sigmax, sigmay, identity
            >>> basis = [ket([0]), ket([1])]
            >>> gate = sigmay()
            >>> gate
            Quantum object: dims = [[2], [2]], shape = (2, 2), type = oper, isherm = True
            Qobj data =
            [[0.+0.j 0.-1.j]
             [0.+1.j 0.+0.j]]
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

            >>> objectives = gate_objectives(
            ...     basis, qutip.gates.cnot(), H, local_invariants=True
            ... )
            >>> for i in range(4):
            ...     assert objectives[i] == Objective(
            ...        initial_state=bell_basis(basis)[i],
            ...        target=qutip.gates.cnot(),
            ...        H=H
            ...     )
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
    target_states = [
        sum([gate[i, j] * basis_states[i] for i in range(gate.shape[0])])
        for j in range(gate.shape[1])
    ]
    return [
        Objective(
            initial_state=initial_state,
            target=target_state,
            H=H,
            c_ops=c_ops,
        )
        for (initial_state, target_state) in zip(basis_states, target_states)
    ]


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


def ensemble_objectives(objectives, Hs):
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

    Returns:
        list[Objective]: List of $n (m+1)$ new objectives that consists of the
        original objectives, plus one copy of the original objectives per
        element of `Hs` where the `H` attribute of each objectives is
        replaced by that element.
    """
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
