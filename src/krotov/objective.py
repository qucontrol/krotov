import sys
import copy
from functools import partial

import attr
import qutip
from qutip.solver import Result as QutipSolverResult
import numpy as np

from .structural_conversions import (
    _nested_list_shallow_copy, extract_controls, extract_controls_mapping,
    _tlist_midpoints, plug_in_pulse_values, control_onto_interval)

__all__ = ['Objective']


#: Workaround for `QuTiP issue 932`_.
#: If True, and only when running on macOS, in :meth:`Objective.mesolve`,
#: replace any array controls with an equivalent function. This results in a
#: signficant slowdown of the propagation, as it circumvents the use of Cython.
#:
#: .. _QuTiP issue 932: https://github.com/qutip/qutip/issues/932
FIX_QUTIP_932 = True


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
    else:
        return op.dag()


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
        # reference.
        return Objective(
            H=_nested_list_shallow_copy(self.H),
            initial_state=self.initial_state,
            target_state=self.target_state,
            c_ops=[_nested_list_shallow_copy(c) for c in self.c_ops])

    @property
    def adjoint(self):
        """The :class:`Objective` containing the adjoint of all components.

        This does not affect the controls in the `H` attribute.
        """
        return Objective(
            H=_adjoint(self.H),
            initial_state=_adjoint(self.initial_state),
            target_state=_adjoint(self.target_state),
            c_ops=[_adjoint(op) for op in self.c_ops])

    def mesolve(self, tlist, rho0=None, e_ops=None, **kwargs):
        """Run :func:`qutip.mesolve.mesolve` on the system of the objective

        Solve the dynamics for the `H` and `c_ops` of the objective. If `rho0`
        is not given, the `initial_state` will be propagated. All other
        arguments will be passed to :func:`qutip.mesolve.mesolve`.

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
                    c_op, controls, mapping[ic+1], tlist)
                for (ic, c_op) in enumerate(self.c_ops)]
        return qutip.mesolve(
            H=H, rho0=rho0, tlist=tlist, c_ops=c_ops, e_ops=e_ops,
            **kwargs)

    def propagate(self, tlist, *, propagator, rho0=None, e_ops=None):
        """Propagates the system of the objective over the entire time grid

        Solve the dynamics for the `H` and `c_ops` of the objective. If `rho0`
        is not given, the `initial_state` will be propagated. This is similar
        to the :meth:`mesolve` method, but instead of using
        :func:`qutip.mesolve.mesolve`, the `propagate` function is used to go
        between time steps. This function is the same as what is passed to
        :func:`.optimize_pulses`. The crucial difference between this and
        :meth:`mesolve` is in the time discretization convention. While
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
        tlist_midpoints = _tlist_midpoints(tlist)
        pulses = [  # defined on the tlist intervals
            control_onto_interval(control, tlist, tlist_midpoints)
            for control in controls]
        for time_index in range(len(tlist)-1):  # index over intervals
            H = plug_in_pulse_values(self.H, pulses, mapping[0], time_index)
            c_ops = [
                plug_in_pulse_values(c_op, pulses, mapping[ic+1], time_index)
                for (ic, c_op) in enumerate(self.c_ops)]
            dt = tlist[time_index+1] - tlist[time_index]
            state = propagator(H, state, dt, c_ops)
            if len(e_ops) == 0:
                result.states.append(state)
            else:
                for (i, oper) in enumerate(e_ops):
                    result.expect[i].append(qutip.expect(oper, state))
        return result


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
