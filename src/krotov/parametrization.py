r"""Classes to realized parametrized optimization pulses."""
import sys
import warnings
from abc import ABCMeta, abstractmethod

import numpy as np


class ParametrizedFunction(metaclass=ABCMeta):
    """Wrapped function, adding a `parametrization` attribute."""

    def __init__(self, func, parametrization):
        self._func = func
        self.parametrization = parametrization

    def __call__(self, t, args):
        return self._func(t, args)


class ParametrizedArray(np.ndarray):
    """Wrapped numpy array, adding a `parametrization` attribute."""

    # See https://numpy.org/doc/stable/user/basics.subclassing.html
    def __new__(cls, input_array, parametrization):
        obj = np.asarray(input_array).view(cls)
        obj.parametrization = parametrization
        if not isinstance(obj.parametrization, Parametrization):
            raise ValueError(
                "parametrization must be a Parametrization instance, not %r"
                % type(parametrization)
            )
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.parametrization = getattr(obj, 'parametrization', None)


class Parametrization(metaclass=ABCMeta):
    """Abstract base class for a parametrizations."""

    @abstractmethod
    def parametrize(self, eps_val):
        return NotImplementedError

    @abstractmethod
    def unparametrize(self, u_val):
        return NotImplementedError

    @abstractmethod
    def derivative(self):
        return NotImplementedError


class TanhParametrization(Parametrization):
    def __init__(self, *, eps_max, eps_min):
        self.eps_max = eps_max
        self.eps_min = eps_min

    def parametrize(self, eps_val):
        ϵ_max = self.eps_max
        ϵ_min = self.eps_min
        ϵ = eps_val
        if ϵ >= ϵ_max or ϵ <= ϵ_min:
            warnings.warn(
                "Pulse value %r out of range (%r, %r) for %s. "
                "Value will be clipped."
                % (ϵ, ϵ_min, ϵ_max, self.__class__.__name__)
            )
        Δ = ϵ_max - ϵ_min
        a = np.clip(
            2 * ϵ / Δ - (ϵ_max + ϵ_min) / Δ,
            -1 + sys.float_info.epsilon,
            1 - sys.float_info.epsilon,
        )
        u = np.arctanh(a)  # -18.4 < u < 18.4
        return u

    def unparametrize(self, u_val):
        ϵ_max = self.eps_max
        ϵ_min = self.eps_min
        u = u_val
        cp = 0.5 * (ϵ_max + ϵ_min)
        cm = 0.5 * (ϵ_max - ϵ_min)
        ϵ = cm * np.tanh(u) + cp
        return ϵ

    def derivative(self, eps_val):
        ϵ_max = self.eps_max
        ϵ_min = self.eps_min
        ϵ = eps_val
        Δ = ϵ_max - ϵ_min
        a = np.clip(
            2 * ϵ / Δ - (ϵ_max + ϵ_min) / Δ,
            -1 + sys.float_info.epsilon,
            1 - sys.float_info.epsilon,
        )
        u = np.arctanh(a)
        return 0.5 * Δ / np.cosh(u) ** 2


class SquareParametrization(Parametrization):
    def parametrize(self, eps_val):
        if eps_val < 0:
            warnings.warn(
                "Pulse value %r < 0 out of range for %s. Clip to 0."
                % (eps_val, self.__class__.__name__)
            )
            eps_val = 0
        return np.sqrt(eps_val)

    def unparametrize(self, u_val):
        return u_val ** 2

    def derivative(self, eps_val):
        u = self.parametrize(eps_val)
        return 2 * u
