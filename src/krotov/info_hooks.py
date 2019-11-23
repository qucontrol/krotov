"""Routines that can be passed as `info_hook` to :func:`.optimize_pulses`"""
import sys
import time

import numpy as np
from uniseg.graphemecluster import grapheme_clusters


__all__ = ['chain', 'print_debug_information', 'print_table']


def _qobj_nbytes(qobj):
    """Return an estimate for the number of bytes of memory required to store a
    Qobj"""
    return (
        qobj.data.data.nbytes
        + qobj.data.indptr.nbytes
        + qobj.data.indices.nbytes
        + sys.getsizeof(qobj)
        + sum(sys.getsizeof(attr) for attr in qobj.__dict__)
    )


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
    forward_states0,
    guess_pulses,
    optimized_pulses,
    g_a_integrals,
    lambda_vals,
    shape_arrays,
    fw_states_T,
    tlist,
    tau_vals,
    start_time,
    stop_time,
    iteration,
    info_vals,
    shared_data,
    out=sys.stdout
):
    r"""Print full debug information about the current Krotov iteration

    This routine is intended to be passed to :func:`.optimize_pulses` as
    `info_hook`, and it exemplifies the full signature of a routine suitable
    for this purpose.

    Keyword Args:
        objectives (list[Objective]): list of the objectives
        adjoint_objectives (list[Objective]): list of the adjoint objectives
        backward_states (list): If available, for each objective, an array-like
            object containing the states of the Krotov backward propagation.
        forward_states: If available (second order only), for each objective,
            an array-like object containing the forward-propagated states under
            the optimized pulses. None otherwise.
        forward_states0: If available (second order only), for each objective,
            an array-like object containing the forward-propagated states under
            the guess pulses. None otherwise.
        guess_pulses (list[numpy.ndarray]): list of guess pulses
        optimized_pulses (list[numpy.ndarray]): list of optimized pulses
        g_a_integrals (numpy.ndarray): array of values
            :math:`\int_0^T g_a(t) \dd t = \int_0^T \frac{\lambda_a}{S(t)}
            \Abs{\Delta \epsilon(t)}^2 \dd t`, for each pulse $\epsilon(t)$.
            The pulse updates $\Delta \epsilon(t)$ are the differences of the
            `optimized_pulses` and the `guess_pulses` (zero in the zeroth
            iteration that only performs a forward-propagation of the guess
            pulses). The quantity $\int g_a(t) \dd t$ is a very useful measure
            of how much the pulse amplitudes changes in each iteration. This
            tells us whether we've chosen good values for $\lambda_a$. Values
            that are too small cause "pulse explosions" which immediately show
            up in $\int_0^T g_a(t) \dd t$.  Also, whether $\int g_a(t) \dd t$
            is increasing or decreasing between iterations gives an indication
            whether the optimization is "speeding up" or "slowing down", and
            thus whether convergence is reached (negligible pulse updates).
        lambda_vals (numpy.ndarray): for each pulse, the value of the
            $\lambda_a$ parameter
        shape_arrays (list[numpy.ndarray]): for each pulse, the array of
            update-shape values $S(t)$
        fw_states_T (list): for each objective, the forward-propagated state
        tlist (numpy.ndarray): array of time grid values on which the states
            are defined
        tau_vals (numpy.ndarray): for each objective, the complex overlap for
            the target state with the forward-propagated state, or None if no
            target state is defined.
        start_time (float): The time at which the iteration started, in epoch
            seconds
        stop_time (float): The time at which the iteration started, in epoch
            seconds
        iteration (int): The current iteration number. For the initial
            propagation of the guess controls, zero.
        shared_data (dict): Dict of data shared between any
            `modify_params_after_iter` and any `info_hook` functions chained
            together via :func:`chain`.
        info_vals (list): List of the return values of the info_hook from
            previous iterations
        out: An open file handle where to write the information. This parameter
            is not part of the `info_hook` interface, and defaults to stdout.
            Use :func:`functools.partial` to pass a different value.

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
    out.write(
        "    ∫gₐ(t)dt: %s\n" % (", ".join(["%.2e" % v for v in g_a_integrals]))
    )
    out.write("    λₐ: %s\n" % (", ".join(["%.2e" % λ for λ in lambda_vals])))
    MB_per_timeslot = sum(_qobj_nbytes(state) for state in fw_states_T) / (
        1024 ** 2
    )
    out.write("    storage (bw, fw, fw0): ")
    if backward_states is None:
        out.write("None, ")
    else:
        n = len(backward_states)
        out.write(
            "[%d * %s(%d)] (%.1f MB), "
            % (
                n,
                backward_states[0].__class__.__name__,
                len(backward_states[0]),
                len(backward_states[0]) * MB_per_timeslot,
            )
        )
    if forward_states is None:
        out.write("None, ")
    else:
        n = len(forward_states)
        out.write(
            "[%d * %s(%d)] (%.1f MB), "
            % (
                n,
                forward_states[0].__class__.__name__,
                len(forward_states[0]),
                len(forward_states[0]) * MB_per_timeslot,
            )
        )
    if forward_states0 is None:
        out.write("None\n")
    else:
        n = len(forward_states0)
        out.write(
            "[%d * %s(%d)] (%.1f MB)\n"
            % (
                n,
                forward_states0[0].__class__.__name__,
                len(forward_states0[0]),
                len(forward_states0[0]) * MB_per_timeslot,
            )
        )
    out.write(
        "    fw_states_T norm: %s\n"
        % (", ".join(["%f" % state.norm() for state in fw_states_T]))
    )
    if not np.any(tau_vals == None):  # noqa
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
    out.flush()


def _grapheme_len(text):
    """Number of graphemes in `text`

    This is the length of the `text` when printed::
        >>> s = 'Â'
        >>> len(s)
        2
        >>> _grapheme_len(s)
        1
    """
    return len(list(grapheme_clusters(text)))


def _rjust(text, width, fillchar=' '):
    """Right-justify text for a total of `width` graphemes

    The `width` is based on graphemes::

        >>> s = 'Â'
        >>> s.rjust(2)
        'Â'
        >>> _rjust(s, 2)
        ' Â'
    """
    len_text = _grapheme_len(text)
    return fillchar * (width - len_text) + text


def print_table(
    J_T,
    show_g_a_int_per_pulse=False,
    J_T_prev=None,
    unicode=True,
    out=sys.stdout,
):
    r"""Print a tabular overview of the functional values in the iteration

    An example output is:

    .. code-block:: console

        iter.        J_T    ∫gₐ(t)dt          J       ΔJ_T         ΔJ  secs
            0   1.00e+00    0.00e+00   1.00e+00        n/a        n/a     0
            1   7.65e-01    2.33e-02   7.88e-01  -2.35e-01  -2.12e-01     1
            2   5.56e-01    2.07e-02   5.77e-01  -2.09e-01  -1.88e-01     1

    The table has the following columns:

    * iteration number
    * value of the final-time functional $J_T$
    * If `show_g_a_int_per_pulse` is True and there is more than one control
      pulse: one column for each pulse containing the value of :math:`\int_0^T
      \frac{\lambda_{a, i}}{S_i(t)} \Abs{\Delta \epsilon_i(t)}^2 \dd t`
    * The value of :math:`\sum_i \int_0^T g_a(\epsilon_i(t)) \dd t =
      \sum_i \int_0^T \frac{\lambda_{a, i}}{S_i(t)}
      \Abs{\Delta \epsilon_i(t)}^2 \dd t`, or just :math:`\int_0^T
      \frac{\lambda_{a}}{S(t)} \Abs{\Delta \epsilon(t)}^2 \dd t` if there
      is only a single control pulse (as in the above example output).
      This value (respectively the individual values with
      `show_g_a_int_per_pulse`) should always be at least three orders of
      magnitude smaller than the pulse fluence :math:`\sum_i\int_0^T
      \Abs{\epsilon_i(t)}^2 \dd t`.  Larger changes in the pulse amplitude may
      be a sign of a "pulse explosion" due to values for $\lambda_{a,i}$ that
      are too small. Changes in :math:`\sum_i \int_0^T \frac{\lambda_{a,
      i}}{S_i(t)}` are often a better indicator of whether the optimization is
      "speeding up"/"slowing down"/reaching convergence than the values of
      $J_T$.
    * The value of the total functional :math:`J = J_T +
      \sum_i \int_0^T g_a(\epsilon_i(t)) \dd t`
    * The change $\Delta J_T$ in the final time functional compared to the
      previous iteration. This should be a negative value, indicating monotonic
      convergence in a minimization of $J_T$.
    * The change $\Delta J$ in the total functional compared to the previous
      iteration. This is evaluated as :math:`\Delta J = \Delta J_T +
      \sum_i \int_0^T g_a(\epsilon_i(t)) \dd t`. Somewhat counter-intuitively,
      $\Delta J$ does not contain a contribution from the $g_a(t)$ of the
      previous iteration. This is because the $\Delta \epsilon_i(T)$ on which
      $g_a(t)$ depends must be evaluated with respect to the same reference
      field (the guess pulse of the *current* iteration), to that
      $\Delta \epsilon_i(T) = 0$ when evaluated with the optimized pulse of the
      previous iteration (i.e., the same guess pulse of the current iteration).
    * The number of seconds in wallclock time spent on the iteration
    * An indicator ``*`` or ``**`` if there is a loss of monotonic convergence
      in $\Delta J_T$ and/or $\Delta J$. Krotov's method mathematically
      guarantees a negative $\Delta J$ in the continuous limit. Assuming there
      are no errors in the time propagation, or in the `chi_constructor` passed
      to :func:`.optimize_pulses`, a loss of monotonic convergence is due to
      the $\lambda_a$ associated with the pulses (via `pulse_options` in
      :func:`.optimize_pulses`) being too small. In practice, we usually don't
      care too much about a loss of monotonic convergence in $\Delta J$, but a
      loss of convergence in $\Delta J_T$ is a serious sign of trouble. It is
      often associated with sharp discontinuous spikes in the optimized pulses,
      or a dramatic increase in the pulse amplitude.

    Args:
        J_T (callable): A function that extract the value of the final time
            functional from the keyword-arguments passed to the `info_hook`.
        show_g_a_int_per_pulse (bool): If True, print a column with the
            value of :math:`\int_0^T g_a(\epsilon_i(t)) \dd t = \int_0^T
            \frac{\lambda_{a, i}}{S_i(t)} \Abs{\Delta \epsilon_i(t)}^2 \dd t`
            for every pulse $\epsilon_i(t)$. Otherwise, only print the sum over
            those integrals for all pulses.
        J_T_prev (None or callable): A function that extract the value of the
            final time functional *from the previous iteration*. If None, use
            the last values from the `info_vals` passed to the `info_hook`.
        unicode (bool): Whether to use unicode symbols for the column headers.
            Some systems have broken monospace fonts in the Jupyter notebook
            that cause the headers not to line up as intended.
        out: An open file handle where to write the table. Defaults to stdout.
    """

    def get_J_T_prev(**kwargs):
        try:
            return kwargs['info_vals'][-1]
        except IndexError:
            return 0

    if J_T_prev is None:
        J_T_prev = get_J_T_prev

    def info_hook(**kwargs):
        iteration = kwargs['iteration']
        n_pulses = len(kwargs['guess_pulses'])
        if iteration == 0:
            out.write(_rjust("iter.", 7))
            out.write(_rjust("J_T", 11))
            if n_pulses > 1:
                if show_g_a_int_per_pulse:
                    for i in range(n_pulses):
                        if unicode:
                            out.write(_rjust("∫gₐ(ϵ%d)dt" % (i + 1), 12))
                        else:
                            out.write(_rjust("g_a_int_%d" % (i + 1), 12))
                if unicode:
                    out.write(_rjust("∑∫gₐ(t)dt", 12))
                else:
                    out.write(_rjust("g_a_int", 12))
            else:
                if unicode:
                    out.write(_rjust("∫gₐ(t)dt", 12))
                else:
                    out.write(_rjust("g_a_int", 12))
            out.write(_rjust("J", 11))
            if unicode:
                out.write(_rjust("ΔJ_T", 11))
            else:
                out.write(_rjust("Delta_J_T", 11))
            if unicode:
                out.write(_rjust("ΔJ", 11))
            else:
                out.write(_rjust("Delta J", 11))
            out.write(_rjust("secs", 6) + "\n")
        J_T_val = J_T(**kwargs)
        Σgₐdt = np.sum(kwargs['g_a_integrals'])
        J = J_T_val + Σgₐdt
        if iteration > 0:
            J_T_prev_val = J_T_prev(**kwargs)
            ΔJ_T = J_T_val - J_T_prev_val
            ΔJ = ΔJ_T + Σgₐdt
        secs = int(kwargs['stop_time'] - kwargs['start_time'])
        out.write("%7d" % iteration)
        out.write(" %10.2e" % J_T_val)
        if n_pulses > 1:
            if show_g_a_int_per_pulse:
                for i in range(n_pulses):
                    out.write(" %11.2e" % kwargs['g_a_integrals'][i])
        out.write(" %11.2e" % Σgₐdt)
        out.write(" %10.2e" % J)
        if iteration == 0:
            out.write(" %10s" % "n/a")
            out.write(" %10s" % "n/a")
        else:
            out.write(" %10.2e" % ΔJ_T)
            out.write(" %10.2e" % ΔJ)
        out.write(" %5d" % secs)
        if iteration > 0:
            if (ΔJ_T > 0) or (ΔJ > 0):
                out.write(" ")
                if ΔJ_T > 0:
                    out.write("*")
                if ΔJ > 0:
                    out.write("*")
        out.write("\n")
        out.flush()
        return J_T_val

    return info_hook


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
