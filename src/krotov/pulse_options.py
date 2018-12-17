import attr

from .shapes import zero_shape, one_shape

__all__ = ['PulseOptions']


def _shape_val_to_callable(val):
    if val == 1:
        return one_shape
    elif val == 0:
        return zero_shape
    else:
        return val


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
            and switch-off. You may also pass the shapes 1 or 0 for a constant
            shape.
    """
    lambda_a = attr.ib()
    shape = attr.ib(default=1, converter=_shape_val_to_callable)

    @shape.validator
    def _check_shape(self, attribute, value):
        if not callable(value):
            raise ValueError("shape must be a callable")
