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
# fmt: off

__version__ = '0.5.0'

__arxiv__ = '1902.11284'

__citation__ = (
    "M. H. Goerz et al., Krotov: A Python implementation of Krotov's method for quantum optimal control, arXiv:%s (2019)"
    % __arxiv__
)

__bibtex__ = r'''
@article{arxiv1902.11284,
    author = {Michael H. Goerz and Daniel Basilewitsch and Fernando Gago-Encinas and Matthias G. Krauss and Karl P. Horn and Daniel M. Reich and Christiane P. Koch},
    title = {Krotov: A {Python} implementation of {Krotov's} method for quantum optimal control},
    year = {2019},
    journal = {arXiv:1902.11284},
}
'''.strip()

# expose submodules for easy import
from . import (
    convergence,
    conversions,
    functionals,
    info_hooks,
    mu,
    objectives,
    parallelization,
    propagators,
    result,
    second_order,
    shapes,
)
# expose primary classes/functions
from .objectives import Objective, ensemble_objectives, gate_objectives
from .optimize import optimize_pulses
from .result import Result


__all__ = [
    'Objective',
    'Result',
    'ensemble_objectives',
    'gate_objectives',
    'optimize_pulses',
]
