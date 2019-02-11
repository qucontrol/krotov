"""Functionals and `chi_constructor` routines"""
import qutip
import numpy as np

from .second_order import _overlap

__all__ = [
    'f_tau',
    'F_re',
    'J_T_re',
    'chis_re',
    'J_T_hs',
    'chis_hs',
    'F_avg',
    'gate',
    'mapped_basis',
]


def f_tau(states_T, objectives, tau_vals=None):
    r"""Average of the complex overlaps with the target states

    That is,

    .. math::

        f_{\tau} = \frac{1}{N} \sum_{i=1}^{N} w_i \tau_i

    where $\tau_i$ are the elements of `tau_vals`, assumed to be

    .. math::

        \tau_i = \Braket{\Psi_i(T)}{\Psi_i^{\tgt}},

    in Hilbert space, or

    .. math::

        \tau_i = \tr\left[\Op{\rho}_i(T)\Op{\rho}_i^{\tgt}\right]

    in Liouville space, where $\ket{\Psi_i}$ or $\Op{\rho}_i$ are the elements
    of `states_T`, and $\ket{\Psi_i^{\tgt}}$ or $\Op{\rho}^{\tgt}$ are the
    target states from the :attr:`~.Objective.target` attribute of the
    objectives. If `tau_vals` are None, they will be calculated internally.

    $N$ is the number of objectives, and $w_i$ is an optional weight for each
    objective. For any objective that has a (custom) `weight` attribute, the
    $w_i$ is taken from that attribute; otherwise, $w_i = 1$. The weights, if
    present, are not automatically normalized, they are assumed to have values
    such that the resulting $f_{\tau}$ lies in the unit circle of the complex
    plane. Usually, this means that the weights should sum to $N$.
    """
    if tau_vals is None:
        tau_vals = [
            _overlap(psi, obj.target)
            for (psi, obj) in zip(states_T, objectives)
        ]
    res = 0j
    for (obj, τ) in zip(objectives, tau_vals):
        if hasattr(obj, 'weight'):
            res += obj.weight * τ
        else:
            res += τ
    return res / len(objectives)


def F_re(states_T, objectives, tau_vals=None):
    r"""Real-part fidelity

    .. math::

        F_{\text{re}} = \Re[f_{\tau}] \quad\in [-1, 1]
    """
    return f_tau(states_T, objectives, tau_vals).real


def J_T_re(states_T, objectives, tau_vals=None):
    r"""Real-part functional

    .. math::

        J_{T,\text{re}} = 1 - F_{\text{re}} \quad\in [0, 2]

    Note:
        If the `states_T` or the target states are mixed states, it is
        preferable to use :func:`J_T_hs`, as $J_{T,\text{re}}$ may take
        negative values for mixed states.
    """
    return 1 - F_re(states_T, objectives, tau_vals)


def chis_re(states_T, objectives, tau_vals):
    r"""States $\ket{\chi}$ for functional $J_{T,\text{re}}$

    .. math::

        \Ket{\chi_i}
        = -\frac{\partial J_{T,\text{re}}}{\partial \bra{\Psi_i(T)}}
        = \frac{1}{2N} w_i \Ket{\Psi^{\tgt}_i}

    with optional weights $w_i$, cf. :func:`f_tau` (default: :math:`w_i=1`). If
    given, the weights should generally sum to $N$.

    Note: `tau_vals` are ignored.
    """
    c = 1.0 / (2 * len(objectives))
    res = []
    for obj in objectives:
        # `obj.target` is assumed to be the "target state" (gate applied to
        # `initial_state`)
        if hasattr(obj, 'weight'):
            res.append(c * obj.weight * obj.target)
        else:
            res.append(c * obj.target)
    return res


def J_T_hs(states_T, objectives, tau_vals=None):
    r"""Hilbert-Schmidt distance measure functional

    .. math::

        J_{T,\text{hs}}
            = \frac{1}{2N} \sum_{i=1}^{N}
                w_i \Norm{\Op{\rho}_i(T) - \Op{\rho}_i^{\tgt}}_{\text{hs}}^2
        \quad \in [0, 2]

    in Liouville space (using the Hilbert-Schmidt norm), or equivalently with
    $\ket{\Psi_i(T)}$ and $\ket{\Psi_i^{tgt}}$ in Hilbert space. The functional
    is evaluated as

    .. math::

        J_{T,\text{hs}}
            = \frac{1}{2N} \sum_{i=1}^{N} w_i \left(
                \Norm{\Op{\rho}_i(T)}_{\text{hs}}^2
                + \Norm{\Op{\rho}^{\tgt}}_{\text{hs}}^2
                - 2 \Re[\tau_i]
            \right)

    where the $\Op{\rho}_i$ are the elements of `states_T`,
    the $\Op{\rho}_i^{\tgt}$ are the target states from the
    :attr:`~.Objective.target` attribute of the objectives,
    and the $\tau_i$ are the elements of `tau_vals` (which
    will be calculated internally if passed as None).

    The $w_i$ are optional weights, cf. :func:`f_tau`. If
    given, the weights should generally sum to $N$.


    Note:
        For pure states (or Hilbert space states), $J_{T,\text{hs}}$ is
        equivalent to $J_{T,\text{re}}$, cf. :func:`J_T_re`. However, the
        backward-propagated states $\chi_i$ obtained from the two functionals
        (:func:`chis_re` and :func:`chis_hs`) are *not* equivalent.
    """
    if tau_vals is None:
        tau_vals = [
            _overlap(psi, obj.target)
            for (psi, obj) in zip(states_T, objectives)
        ]
    res = 0.0
    hs = 'l2'  # qutip's name for HS-norm for state vectors
    if states_T[0].type == 'oper':
        hs = 'fro'  # qutip's name for HS-norm for density matrices
    for (obj, ρ, τ) in zip(objectives, states_T, tau_vals):
        ρ_tgt = obj.target
        if hasattr(obj, 'weight'):
            res += obj.weight * (
                ρ.norm(hs) ** 2 + ρ_tgt.norm(hs) ** 2 - 2 * τ.real
            )
        else:
            res += ρ.norm(hs) ** 2 + ρ_tgt.norm(hs) ** 2 - 2 * τ.real
    return res / (2 * len(objectives))


def chis_hs(states_T, objectives, tau_vals):
    r"""States $\Op{\chi}$ for functional $J_{T,\text{hs}}$

    .. math::

        \Op{\chi}_i
        = -\frac{\partial J_{T,\text{sm}}}
                {\partial \langle\!\langle \Op{\rho}_i(T)\vert}
        = \frac{1}{2N} w_i
          \left(\Op{\rho}^{\tgt}_i - \Op{\rho}_i(T)\right)

    with optional weights $w_i$, cf. :func:`f_tau` (default: :math:`w_i=1`).

    This is derived from $J_{T,\text{hs}}$ rewritten in the abstract
    Hilbert-Schmidt notation :math:`\langle\!\langle a \vert b \rangle\!\rangle
    \equiv \tr[a^\dagger b]`:

    .. math::

        J_{T,\text{hs}} = \frac{-1}{2N} \sum_{i=1}^{N}  w_i \big(
            \underbrace{
            \langle\!\langle \Op{\rho}_i(T) \vert
                \Op{\rho}_i^{\tgt} \rangle\!\rangle
            + \langle\!\langle \Op{\rho}_i^{\tgt}\vert
                 \Op{\rho}_i(T) \rangle\!\rangle
            }_{=2\Re[\tau_i]}
            - \underbrace{
              \langle\!\langle \Op{\rho}_i(T) \vert
                \Op{\rho}_i(T) \rangle\!\rangle
            }_{=\Norm{\Op{\rho}_i(T)}_{\text{hs}}^2}
            - \underbrace{
              \langle\!\langle \Op{\rho}_i^{\tgt} \vert
                \Op{\rho}_i^{\tgt} \rangle\!\rangle
            }_{=\Norm{\Op{\rho}^{\tgt}}_{\text{hs}}^2}
        \big).

    Note that `tau_vals` is ignored.
    """
    c = 1.0 / (2 * len(objectives))
    res = []
    for (obj, ρ) in zip(objectives, states_T):
        ρ_tgt = obj.target
        if hasattr(obj, 'weight'):
            w = obj.weight
            res.append(c * w * (ρ_tgt - ρ))
        else:
            res.append(c * (ρ_tgt - ρ))
    return res


def F_avg(states_T, basis_states, gate, mapped_basis_states=None):
    r"""Average gate fidelity

    .. math::

        F_{\text{avg}}
            = \int \big\langle \Psi \big\vert
                \Op{O}^\dagger \DynMap[\ketbra{\Psi}{\Psi}] \Op{O} \big\vert
                \Psi \big\rangle \dd \Psi

    where $\Op{O}$ is the target `gate`, and $\DynMap$ represents the dynamical
    map from time zero to $T$.

    In Liouville space, this is numerically evaluated as

    .. math::

        F_{\text{avg}} = \frac{1}{N (N+1)}
            \sum_{i,j=1}^N \left(
                \big\langle
                    \phi_i \big\vert
                    \Op{O}^\dagger \Op{\rho}_{ij} \Op{O} \big\vert
                    \phi_j
                \big\rangle
                + \big\langle
                    \phi_i \big\vert
                    \Op{O}^\dagger  \Op{\rho}_{jj} \Op{O} \big\vert
                    \phi_i
                \big\rangle
            \right),

    where :math:`\ket{\phi_i}` is the :math:`i`'th element of `basis_states`,
    and :math:`\Op{\rho}_{ij}` is the :math:`(i-1) N + j`'th element of
    `states_T`, that is, :math:`\Op{\rho}_{ij} =
    \DynMap[\ketbra{\phi_i}{\phi_j}]`, with :math:`N` the dimension of the
    Hilbert space.

    In Hilbert space (unitary dynamics), this simplifies to

    .. math::

        F_{\text{avg}} = \frac{1}{N (N+1)} \left(
                \Abs{\tr\left[\Op{O}^\dagger \Op{U}\right]}^2
                + \tr\left[\Op{O}^\dagger \Op{U} \Op{U}^\dagger \Op{O}\right]
            \right),

    where $\Op{U}$ the gate that maps `basis_states` to the result of a forward
    propagation of those basis states, stored in `states_T`.

    Args:
        states_T (list[qutip.Qobj]): The forward propagated states. For
            dissipative dynamics, this must be the forward propagation of the
            full basis of Liouville space, that is, all $N^2$ dyadic
            combinations of the Hilbert space logical basis states.
            For unitary dynamics, the $N$ forward-propagated `basis_states`.
        basis_states (list[qutip.Qobj]): The $N$ Hilbert space logical basis
            states
        gate (qutip.Qobj): The $N \times N$ quantum gate in the logical
            subspace, e.g. :func:`qutip.qip.gates.cnot()`.
        mapped_basis_states (None or list[qutip.Qobj]): If given, the result of
            applying gate to `basis_states`. If not given, this will be
            calculated internally via :func:`mapped_basis`. It is recommended
            to pass pre-calculated `mapped_basis_states` when evaluating
            $F_{\text{avg}}$ repeatedly for the same target.
    """
    N = len(basis_states)
    if gate.shape != (N, N):
        raise ValueError(
            "Shape of gate is incompatible with number of basis states"
        )
    if states_T[0].type == 'oper':
        if len(states_T) != N * N:
            raise ValueError(
                "Evaluating F_avg for density matrices requires %d states "
                "(forward-propagation of all dyadic combinations of "
                "%d basis states), not %d" % (N * N, N, len(states_T))
            )
        return _F_avg_rho(states_T, basis_states, gate, mapped_basis_states)
    elif states_T[0].type == 'ket':
        if len(states_T) != N:
            raise ValueError(
                "Evaluating F_avg for hilbert space states requires %d states "
                "(forward-propagation of all basis states), not %d"
                % (N, len(states_T))
            )
        return _F_avg_psi(states_T, basis_states, gate)
    else:
        raise ValueError("Invalid type of state: %s" % states_T[0].type)


def _F_avg_rho(states_T, basis_states, gate, mapped_basis_states):
    """Implementation of F_avg in Liouville space"""
    if mapped_basis_states is None:
        mapped_basis_states = mapped_basis(gate, basis_states)
    N = len(basis_states)
    F = 0
    for j in range(N):
        ρ_jj = states_T[j * N + j]
        Oϕ_j = mapped_basis_states[j]
        for i in range(N):
            ρ_ij = states_T[i * N + j]
            Oϕ_i = mapped_basis_states[i]
            F += _overlap(Oϕ_i, ρ_ij(Oϕ_j)) + _overlap(Oϕ_i, ρ_jj(Oϕ_i))
    assert abs(F.imag) < 1e-14
    return F.real / (N * (N + 1))


def _F_avg_psi(states_T, basis_states, O):
    """Implementation of F_avg in Hilbert space"""
    N = len(basis_states)
    U = gate(basis_states, states_T)
    F = abs((O.dag() * U).tr()) ** 2 + (O.dag() * U * U.dag() * O).tr()
    assert abs(F.imag) < 1e-14
    return F.real / (N * (N + 1))


def gate(basis_states, states_T):
    """Gate that maps `basis_states` to `states_T`

    Example:

        >>> from qutip import ket
        >>> basis = [ket(nums) for nums in [(0, 0), (0, 1), (1, 0), (1, 1)]]
        >>> states_T = mapped_basis(qutip.gates.cnot(), basis)
        >>> U = gate(basis, states_T)
        >>> assert (U - qutip.gates.cnot()).norm() < 1e-15
    """
    N = len(basis_states)
    U = np.zeros((N, N), dtype=np.complex128)
    for j in range(N):
        for i in range(N):
            U[i, j] = basis_states[i].overlap(states_T[j])
    dims = [basis_states[0].dims[0], states_T[0].dims[0]]
    return qutip.Qobj(U, dims=dims)


def mapped_basis(gate, basis_states):
    """Result of applying `gate` to `basis_states`

    Example:

        >>> from qutip import ket
        >>> basis = [ket(nums) for nums in [(0, 0), (0, 1), (1, 0), (1, 1)]]
        >>> states = mapped_basis(qutip.gates.cnot(), basis)
        >>> assert (states[0] - ket((0,0))).norm() < 1e-15
        >>> assert (states[1] - ket((0,1))).norm() < 1e-15
        >>> assert (states[2] - ket((1,1))).norm() < 1e-15  # swap (1, 1) ...
        >>> assert (states[3] - ket((1,0))).norm() < 1e-15  # ... and (1, 0)
    """
    return tuple(
        [
            sum([gate[i, j] * basis_states[i] for i in range(gate.shape[0])])
            for j in range(gate.shape[1])
        ]
    )
