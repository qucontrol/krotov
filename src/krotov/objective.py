import copy

import attr

__all__ = ['Objective']


@attr.s
class Objective():
    """A single objective for optimization with Krotov's method

    Attributes:
        initial_state (qutip.Qobj): The initial state
        target_state (qutip.Qobj): The desired target state
        H (list): The time-dependent Hamiltonian,
            cf. :func:`qutip.mesolve.mesolve`. This includes the control
            fields.
        c_ops (list):  List of collapse operators,
            cf. :func:`~qutip.mesolve.mesolve`.
    """
    H = attr.ib()
    initial_state = attr.ib()
    target_state = attr.ib()
    c_ops = attr.ib(default=[])

    def __copy__(self):
        # When we use copy.copy(objective), we want a
        # semi-deep copy where nested lists in the Hamiltonian and the c_ops
        # are re-created (copy by value), but non-list elements are copied by
        # reference. The _extract_controls relies heavily on this exact
        # behavior.
        return Objective(
            H=_nested_list_shallow_copy(self.H),
            initial_state=self.initial_state,
            target_state=self.target_state,
            c_ops=[_nested_list_shallow_copy(c) for c in self.c_ops])


def _nested_list_shallow_copy(l):
    if isinstance(l, list):
        return [copy.copy(h) if isinstance(h, list) else h for h in l]
    else:
        return l
