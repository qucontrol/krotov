import attr

__all__ = ['PulseOptions']


@attr.s
class PulseOptions():
    """Options for the optimization of a control pulse

    Attributes:
        lambda_a (float): Krotov step size. This governs the overall magnitude
            of the pulse update. Large values result in small updates. Small
            values may lead to sharp spikes and numerical instability.
        shape (callable): Function S(t) in the range [0, 1] that scales the
            pulse update for the pulse value at t. This can be used to ensure
            boundary conditions (S(0) = S(T) = 0), and enforce smooth switch-on
            and switch-off
        filter (callable or None): A function that manipulates the pulse after
            each OCT iteration, e.g. by applying a spectral filter.
    """
    lambda_a = attr.ib()
    shape = attr.ib(default=lambda t: 1)
    filter = attr.ib(default=None)
