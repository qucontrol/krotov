"""Routines that can be passed as `info_hook` to :func:`.optimize_pulses`"""
import sys
import time
import numpy as np

__all__ = ['chain', 'print_debug_information']


def chain(*hooks):
    """Chain multiple `info_hook` callables together"""

    def info_hook(**kwargs):
        result = []
        for hook in hooks:
            res = hook(**kwargs)
            if res is not None:
                result.append(result)
        if len(result) > 0:
            return tuple(result)
        else:
            return None

    return info_hook


def print_debug_information(
        objectives, adjoint_objectives, backward_states, forward_states,
        optimized_pulses, lambda_vals, shape_arrays, fw_states_T, tau_vals,
        start_time, stop_time, iteration, out=sys.stdout):
    """Print full debug information about the current Krotov iteration

    This routine is intended to be passed to :func:`.optimize_pulses` as
    `info_hook`, and it exemplifies the full signature of a routine suitable
    for this purpose.

    Args:
        objectives (list[Objective]): list of the objectives
        adjoint_objectives (list[Objective]): list of the adjoint objectives
        backward_states (list): If available, for each objective, an array-like
            object containing the states of the Krotov backward propagation.
        forward_states: If available, for each objective, an array-like object
            containing the forward-propagated states under the optimized
            pulses.
        optimized_pulses (list): list of optimized pulses, where each pulse is
            a numpy array
        lambda_vals (list): for each pulse, the value of the $\lambda_a$
            parameter
        shape_arrays (list): for each pulse, the array of update-shape values
            $S(t)$
        fw_states_T (list): for each objective, the forward-propagated state
        tau_vals (list): for each objective, the complex overlap for the
            forward-propagated state with the target state
        start_time (float): The time at which the iteration started, in epoch
            seconds
        stop_time (float): The time at which the iteration started, in epoch
            seconds
        iteration (int): The current iteration number. For the initial
            propagation of the guess controls, 0.
        out: An open file handle where to write the information.

    Note:
        This routine implements the full signature of an `info_hook` in
        :func:`.optimize_pulses`. However, since :func:`.optimize_pulses` will
        always call the `info_hook` with keyword arguments, it is usually much
        simpler to define info-hooks using Python's variable keyword arguments
        syntax (``**args``). For example, consider the following `info_hook`
        that prints (and stores) the value of the real-part-gate-fidelity::

            def print_fidelity(**args):
                F_re = np.average(np.array(args['tau_vals']).real)
                print("    F = %f" % F_re)
                return F_re
    """
    out.write("Iteration %d\n" % iteration)
    if iteration == 0:
        out.write("    objectives:\n")
        for (i, obj) in enumerate(objectives):
            out.write("        %d:%s\n" % (i+1, obj))
        out.write("    adjoint objectives:\n")
        for (i, obj) in enumerate(adjoint_objectives):
            out.write("        %d:%s\n" % (i+1, obj))
        out.write("    λₐ: %s\n" % (
            ", ".join(["%.2e" % λ for λ in lambda_vals])))
    out.write("    duration: %.1f secs (started at %s)\n" % (
        (stop_time-start_time),
        time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))))
    out.write("    τ: %s\n" % (", ".join([
        "(%.2e:%.2fπ)" % (abs(z), np.angle(z)/np.pi) for z in tau_vals])))
