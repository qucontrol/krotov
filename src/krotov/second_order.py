"""Support functions for the second-order update equation"""
from abc import ABC, abstractmethod
from typing import Optional


__all__ = ['Sigma', 'numerical_estimate_A']


class Sigma(ABC):
    """Function σ(t) for the second order update equation.

    This is an abstract bases class. For any optimization that requires the
    second-order update equation, an appropriate problem-specific subclass of
    :class:`Sigma` must be implemented that defines

    * the evaluation of σ(t) in :meth:`__call__`

    * the update of any values that σ(t) depends on parametrically (typically:
      any of the parameters A, B, C), in :meth:`refresh`.

    An instantiation of that subclass is then passed as `sigma` to
    :func:`.optimize_pulses`.
    """

    @abstractmethod
    def __call__(self, t):  # pragma: nocover
        """Evaluate σ(t)"""
        raise NotImplementedError()

    @abstractmethod
    def refresh(
        self,
        forward_states,
        forward_states0,
        chi_states,
        chi_norms,
        optimized_pulses,
        guess_pulses,
        objectives,
        result,
    ):  # pragma: nocover
        """Recalculate the parametric dependencies of σ(t)

        This is called at the end of each control iteration, and may be used to
        estimate the internal parameters in σ(t)

        Args:
            forward_states (list): For each objective, an array-like container
                (cf. `storage` in :func:`.optimize_pulses`) of the initial
                state forward-propagated under optimized controls from the
                current iteration.
            forward_states0 (list): The forward-propagated states under the
                guess controls of the current iteration.
            chi_states (list): The (normalized) boundary condition for the
                backward-propagation in the current iteration, as returned by
                the `chi_constructor` argument to :func:`.optimize_pulses`.
            chi_norms (list): The norms of the un-normalized `chi_states`.
            optimized_pulses (list[numpy.ndarray]) list of optimized pulses
                from the current iteration
            guess_pulses (list[numpy.ndarray]) list of guess pulses for the
                current iteration
            objectives (list[Objective]): The control objectives
            result (Result): The result object, up-to-date for the current
                iteration
        """
        raise NotImplementedError()


def _overlap(a, b) -> Optional[complex]:
    """Complex overlap of two quantum objects.

    If `a`, `b` are not quantum objects or are not compatible, return None.
    """
    try:
        if a.type == b.type == 'oper':
            if a.isherm:
                return complex((a * b).tr())
            else:
                return complex((a.dag() * b).tr())
        else:
            return a.overlap(b)
    except AttributeError:
        return None


def numerical_estimate_A(
    forward_states, forward_states0, chi_states, chi_norms, Delta_J_T
):
    r"""Update the second-order parameter $A$.

    Calculate the new value of $A$ according to the equation

    .. math::

        A^{(i+1)} = \frac{
            \sum_k 2 \Re \Braket{\chi_k(T)}{\Delta\phi_k(T)} + \Delta J_T
        }{
            \sum_k \Braket{\Delta \phi_k(T)}{\Delta\phi_k(T)}
        },

    where $\Delta\phi_k$ is the difference of the `forward_states`
    $\ket{\phi_k^{(i)}}$ propagated under the optimized pulse of iteration
    $(i)$, and the `forward_states0` $\ket{\phi_k^{(i-1)}}$ propagated under
    the guess pulse of iteration $(i)$ -- that is, the guess pulse of iteration
    $(i-1)$; and $\Delta J_T$ is the difference of the final time functional,

    .. math::

        \Delta J_T
        = J_T(\{\ket{\phi_k^{(i)}(T)}\} - J_T(\{\ket{\phi_k^{(i-1)}(T)}\}.

    Args:
        forward_states (list): For each objective, the result of a
            forward-propagation with the optimized pulses of the current
            iteration.
        forward_states0 (list): For each objective, the result of a
            forward-propagation with the guess pulses of the current iteration
        chi_states (list): For each objective, the normalized boundary state
            $\ket{\chi_k(T)}/\Abs{\ket{\chi_k(T)}}$ for the
            backward-propagation with the guess pulse of the current iteration.
        chi_norms (list): The norms of the `chi_states`
        Delta_J_T (float): The value by which the final time functional
            improved in the current iteration.
    """
    n = len(forward_states0)  # the number of objectives
    Δϕ = [forward_states[k][-1] - forward_states0[k][-1] for k in range(n)]
    Δϕ_nrmsq = [_overlap(Δϕ[k], Δϕ[k]).real for k in range(n)]
    denom = sum(Δϕ_nrmsq)
    if denom > 1.0e-30:
        numer = (
            sum(
                [
                    (2 * chi_norms[k] * _overlap(chi_states[k], Δϕ[k])).real
                    for k in range(n)
                ]
            )
            + Delta_J_T
        )
        return numer / denom
    else:
        return 0
