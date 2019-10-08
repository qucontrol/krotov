import os
import qutip
import krotov
import mpl


def plot_iterations(opt_result):
    """Plot the control fields in population dynamics over all iterations.

    This depends on ``store_all_pulses=True`` in the call to
    `optimize_pulses`.
    """
    fig_width = 15  # cm
    fig_height = 4.5  # cm
    left_margin = 1.25  # cm
    right_margin = 0.15  # cm
    top_margin = 0.15  # cm
    bottom_margin = 0.75  # cm
    gap = 1.25  # cm

    w = 0.5 * (fig_width - left_margin - right_margin - gap)
    h = fig_height - bottom_margin - top_margin

    fig = mpl.new_figure(fig_width, fig_height)
    ax_ctr = fig.add_axes(
        [
            left_margin / fig_width,
            bottom_margin / fig_height,
            w / fig_width,
            h / fig_height,
        ]
    )
    ax_dyn = fig.add_axes(
        [
            (left_margin + w + gap) / fig_width,
            bottom_margin / fig_height,
            w / fig_width,
            h / fig_height,
        ]
    )

    proj0 = qutip.ket2dm(qutip.ket("0"))
    proj1 = qutip.ket2dm(qutip.ket("1"))

    n_iters = len(opt_result.iters)
    for (iteration, pulses) in zip(opt_result.iters, opt_result.all_pulses):
        controls = [
            krotov.conversions.pulse_onto_tlist(pulse)
            for pulse in pulses
        ]
        objectives = opt_result.objectives_with_controls(controls)
        dynamics = objectives[0].mesolve(
            opt_result.tlist, e_ops=[proj0, proj1]
        )
        if iteration == 0:
            ls = '--'  # dashed
            alpha = 1  # full opacity
            ctr_label = 'guess'
            pop_labels = ['0 (guess)', '1 (guess)']
        elif iteration == opt_result.iters[-1]:
            ls = '-'  # solid
            alpha = 1  # full opacity
            ctr_label = 'optimized'
            pop_labels = ['0 (opt.)', '1 (opt.)']
        else:
            ls = '-'  # solid
            alpha = 0.5 * float(iteration) / float(n_iters)  # max 50%
            ctr_label = None
            pop_labels = [None, None]
        ax_ctr.plot(
            dynamics.times,
            controls[0],
            label=ctr_label,
            color='black',
            ls=ls,
            alpha=alpha,
        )
        ax_dyn.plot(
            dynamics.times,
            dynamics.expect[0],
            label=pop_labels[0],
            color='#1f77b4',  # default blue
            ls=ls,
            alpha=alpha,
        )
        ax_dyn.plot(
            dynamics.times,
            dynamics.expect[1],
            label=pop_labels[1],
            color='#ff7f0e',  # default orange
            ls=ls,
            alpha=alpha,
        )
    ax_ctr.legend()
    ax_ctr.set_xlabel('time')
    ax_ctr.set_ylabel('control amplitude')
    ax_dyn.legend()
    ax_dyn.set_xlabel('time')
    ax_dyn.set_ylabel('population')
    mpl.add_panel_labels([ax_ctr, ax_dyn], fmt='({label})', offset_y=4)
    return fig


def main():
    scriptfile = os.path.abspath(__file__)
    figdir = os.path.dirname(scriptfile)
    dumpfile = os.path.join(figdir, '../examples/tls_state_to_state.dump')
    outfile = os.path.splitext(scriptfile)[0] + '.pdf'
    opt_result = krotov.result.Result.load(dumpfile)
    fig = plot_iterations(opt_result)
    fig.savefig(outfile, transparent=True)


if __name__ == "__main__":
    main()
