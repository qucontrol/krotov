"""The main function exposed here is :func:`.optimize_pulses`.

This acts on a list of :class:`.Objective` instances, which may either be
constructed manually, or through the helper functions :func:`.gate_objectives`
and :func:`.ensemble_objectives`.

The submodules contain various auxiliary functions for constructing arguments
for :func:`.optimize_pulses`, for common use cases. This includes
functions for constructing control fields (guess pulses), optimization
functionals, convergence checks, analysis tools, and propagators, as well as
more technical routines for parallelization, low-level data conversion, and
estimators for second-order updates.
"""

__version__ = '0.2.0+dev'

# expose submodules for easy import
from . import shapes
from . import structural_conversions
from . import propagators
from . import functionals
from . import info_hooks
from . import objectives
from . import mu
from . import result
from . import convergence
from . import second_order
from . import parallelization

# expose primary classes/functions
from .objectives import Objective, gate_objectives, ensemble_objectives
from .optimize import optimize_pulses

__all__ = [
    'optimize_pulses',
    'Objective',
    'gate_objectives',
    'ensemble_objectives',
]
