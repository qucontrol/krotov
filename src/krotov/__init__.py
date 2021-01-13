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

__version__ = '1.2.1'

__arxiv__ = '1902.11284'

__citation__ = (
    "M. H. Goerz et al., "
    "Krotov: A Python implementation of Krotov's method for quantum optimal "
    "control, "
    "SciPost Phys. 7, 080 (2019)"
)

__bibtex__ = r'''
@article{GoerzSPP2019,
    author = {Michael H. Goerz and Daniel Basilewitsch and Fernando Gago-Encinas and Matthias G. Krauss and Karl P. Horn and Daniel M. Reich and Christiane P. Koch},
    title = {Krotov: A {Python} implementation of {Krotov's} method for quantum optimal control},
    journal={SciPost Phys.},
    volume={7},
    pages={80},
    year={2019},
    doi={10.21468/SciPostPhys.7.6.080},
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
