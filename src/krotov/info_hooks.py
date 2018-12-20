"""Routines that can be passed as `info_hook` to :func:`.optimize_pulses`"""
import sys
import time
import numpy as np

__all__ = ['chain', 'print_debug_information']


def chain(*hooks):
    """Chain multiple `info_hook` or `modify_params_after_iter` callables
    together.

    Example:

        >>> def print_fidelity(**kwargs):
        ...     F_re = np.average(np.array(kwargs['tau_vals']).real)
        ...     print("    F = %f" % F_re)
        >>> info_hook = chain(print_debug_information, print_fidelity)

    Note:

        Functions that are connected via :func:`chain` may use the
        `shared_data` share the same `shared_data` argument, which they can use
        to communicate down the chain.
    """

    def info_hook(**kwargs):
        result = []
        for hook in hooks:
            res = hook(**kwargs)
            if res is not None:
                result.append(res)
        if len(result) > 0:
            if len(result) == 1:
                return result[0]
            else:
                return tuple(result)
        else:
            return None

    return info_hook


def print_debug_information(
    *,
    objectives,
    adjoint_objectives,
    backward_states,
    forward_states,
    optimized_pulses,
    lambda_vals,
    shape_arrays,
    fw_states_T,
    tau_vals,
    start_time,
    stop_time,
    iteration,
    shared_data,
    out=sys.stdout
):
    """Print full debug information about the current Krotov iteration

    This routine is intended to be passed to :func:`.optimize_pulses` as
    `info_hook`, and it exemplifies the full signature of a routine suitable
    for this purpose.

    Keyword Args:
        objectives (list[Objective]): list of the objectives
        adjoint_objectives (list[Objective]): list of the adjoint objectives
        backward_states (list): If available, for each objective, an array-like
            object containing the states of the Krotov backward propagation.
        forward_states: If available, for each objective, an array-like object
            containing the forward-propagated states under the optimized
            pulses.
        optimized_pulses (list[numpy.ndarray]): list of optimized pulses
        lambda_vals (list[float]): for each pulse, the value of the $\lambda_a$
            parameter
        shape_arrays (list[numpy.ndarray]): for each pulse, the array of
            update-shape values $S(t)$
        fw_states_T (list): for each objective, the forward-propagated state
        tau_vals (list[complex]): for each objective, the complex overlap for
            the forward-propagated state with the target state
        start_time (float): The time at which the iteration started, in epoch
            seconds
        stop_time (float): The time at which the iteration started, in epoch
            seconds
        iteration (int): The current iteration number. For the initial
            propagation of the guess controls, 0.
        shared_data (dict): Dict of data shared between any
            `modify_params_after_iter` and any `info_hook` functions chained
            together via :func:`chain`.
        out: An open file handle where to write the information. This parameter
            is not part of the `info_hook` interface. Use
            :func:`functools.partial` to pass a different value.

    Note:
        This routine implements the full signature of an `info_hook` in
        :func:`.optimize_pulses`, excluding `out`. However, since
        the `info_hook` only allows for keyword arguments, it is usually much
        simpler to use Python's variable keyword arguments syntax
        (``**kwargs``). For example, consider the following `info_hook` that
        prints (and stores) the value of the real-part gate fidelity::

            def print_fidelity(**kwargs):
                F_re = np.average(np.array(kwargs['tau_vals']).real)
                print("    F = %f" % F_re)
                return F_re
    """
    out.write("Iteration %d\n" % iteration)
    if iteration == 0:
        out.write("    objectives:\n")
        for (i, obj) in enumerate(objectives):
            out.write("        %d:%s\n" % (i + 1, obj))
        out.write("    adjoint objectives:\n")
        for (i, obj) in enumerate(adjoint_objectives):
            out.write("        %d:%s\n" % (i + 1, obj))
        out.write(
            "    λₐ: %s\n" % (", ".join(["%.2e" % λ for λ in lambda_vals]))
        )
        out.write(
            "    S(t) (ranges): %s\n"
            % (
                ", ".join(
                    ["[%f, %f]" % (np.min(S), np.max(S)) for S in shape_arrays]
                )
            )
        )
    out.write(
        "    duration: %.1f secs (started at %s)\n"
        % (
            (stop_time - start_time),
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time)),
        )
    )
    out.write(
        "    optimized pulses (ranges): %s\n"
        % (", ".join([_pulse_range(pulse) for pulse in optimized_pulses]))
    )
    if backward_states is None:
        out.write("    backward states: None\n")
    else:
        n = len(backward_states)
        out.write(
            "    backward states: [%d * %s(%d)]\n"
            % (
                n,
                backward_states[0].__class__.__name__,
                len(backward_states[0]),
            )
        )
    if forward_states is None:
        out.write("    forward states: None\n")
    else:
        n = len(forward_states)
        out.write(
            "    forward states: [%d * %s(%d)]\n"
            % (n, forward_states[0].__class__.__name__, len(forward_states[0]))
        )
    out.write(
        "    fw_states_T norm: %s\n"
        % (", ".join(["%f" % state.norm() for state in fw_states_T]))
    )
    out.write(
        "    τ: %s\n"
        % (
            ", ".join(
                [
                    "(%.2e:%.2fπ)" % (abs(z), np.angle(z) / np.pi)
                    for z in tau_vals
                ]
            )
        )
    )


def _pulse_range(pulse):
    """Return the range information about the given pulse value array

    Example:

        >>> _pulse_range(np.array([-1, 1, 5]))
        '[-1.00, 5.00]'
        >>> _pulse_range(np.array([-1, 1+5j, 5]))
        '[(r:-1.00, i:0.00), (r:5.00, i:5.00)]'
    """
    if np.iscomplexobj(pulse):
        pulse_real = pulse.real
        pulse_imag = pulse.imag
        r0, r1 = np.min(pulse_real), np.max(pulse_real)
        i0, i1 = np.min(pulse_imag), np.max(pulse_imag)
        return "[(r:%.2f, i:%.2f), (r:%.2f, i:%.2f)]" % (r0, i0, r1, i1)
    else:
        return "[%.2f, %.2f]" % (np.min(pulse), np.max(pulse))
