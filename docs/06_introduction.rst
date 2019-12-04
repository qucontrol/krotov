Introduction
============

Quantum information science has changed our perception of quantum
physics from passive understanding to a source of technological
advances :cite:`AcinNJP18`. By way of actively exploiting
the two essential elements of quantum physics, coherence and
entanglement, technologies such as quantum
computing :cite:`NielsenChuang` or quantum
sensing :cite:`DegenRMP17` hold the promise for solving
computationally hard problems or reaching unprecedented sensitivity.
These technologies rely on the ability to accurately perform quantum
operations for increasingly complex quantum systems. Quantum optimal
control allows to address this challenge by providing a set of tools to
devise and implement shapes of external fields that accomplish a given
task in the best way possible :cite:`GlaserEPJD2015`.
Originally developed in the context of molecular
physics :cite:`Tannor92,GrossJCP92` and nuclear magnetic
resonance :cite:`MurdochJMR87,GlaserCPL89`, quantum optimal
control theory has been adapted to the specific needs of quantum
information science in recent
years :cite:`GlaserEPJD2015,KochJPCM16`. Calculation of
optimized external field shapes for tasks such as state preparation or
quantum gate implementation have thus become
standard :cite:`GlaserEPJD2015`, even for large Hilbert
space dimensions as encountered in e.g. Rydberg
atoms :cite:`CuiQST17,PatschPRA18`. Experimental
implementation of the calculated field shapes, using arbitrary waveform
generators, has been eased by the latter becoming available
commercially. Successful demonstration of quantum operations in various
experiments :cite:`GlaserEPJD2015,LovecchioPRA16,vanFrankSciRep16,OfekNat16,SorensenNat16,HeeresNatComm17,HeckPNAS18,FengPRA18,OmranS2019,Larrouy`
attests to the level of maturity that quantum optimal control in quantum
technologies has reached.

In order to calculate optimized external field shapes, two choices need to be
made – about the optimization functional and about the optimization method.
The functional consists of the desired figure of merit, such as a gate or state
preparation error, as well as additional constraints, such as amplitude or
bandwidth restrictions :cite:`GlaserEPJD2015,KochJPCM16`.
Optimal control methods in general can be classified into gradient-free and
gradient-based algorithms that either evaluate the optimization functional
alone or together with its gradient :cite:`GlaserEPJD2015`.
Gradient-based methods typically converge faster, unless the number of
optimization parameters can be kept small.
Most gradient-based methods rely on the iterative solution of a set of coupled
equations that include forward propagation of initial states, backward
propagation of adjoint states, and the control update :cite:`GlaserEPJD2015`.
A popular representative of concurrent update methods is GRadient Ascent Pulse
Engineering (GRAPE) :cite:`KhanejaJMR05`.
Krotov's method, in contrast, requires sequential
updates :cite:`Tannor92,ReichJCP12`.
This comes with the advantage of guaranteed monotonic convergence and
obviates the need for a line search in the direction of the
gradient :cite:`EitanPRA11`.


The choice of Python as an implementation language is due to Python's
easy-to-learn syntax, expressiveness, and immense popularity in the
scientific community. Moreover, the `QuTiP library`_ exists,
providing a general purpose tool to numerically describe quantum systems
and their dynamics. QuTiP already includes basic versions of other
popular quantum control algorithms such as GRAPE and the gradient-free
CRAB :cite:`CanevaPRA11`. The `Jupyter notebook`_ framework is available to provide an ideal
platform for the interactive exploration of the :mod:`krotov` package's
capabilities, and to facilitate reproducible research workflows.

The :mod:`krotov` package targets both students wishing to
enter the field of quantum optimal control, and researchers in the
field. By providing a comprehensive set of :ref:`krotov-example-notebooks`, we
enable users of
our package to explore the formulation of typical control problems, and
to understand how Krotov's method can solve them. These examples are
inspired by recent
publications :cite:`MullerQIP11,GoerzPRA2014,GoerzNJP2014,WattsPRA2015,GoerzPRA2015,BasilewitschNJP2017`,
and thus show the use of the method in the purview of current research.
In particular, the package is not restricted to closed quantum systems,
but can fully address open system dynamics, and thus aide in the
development of Noisy Intermediate-Scale Quantum (NISQ)
technology :cite:`PreskillQ2018`. Optimal control is also
increasingly important in the design of
experiments :cite:`GlaserEPJD2015,LovecchioPRA16,vanFrankSciRep16,OfekNat16,SorensenNat16,HeeresNatComm17,HeckPNAS18,FengPRA18,OmranS2019,Larrouy`,
and we hope that the availability of an easy-to-use implementation of
Krotov's method will facilitate this further.

Large Hilbert space
dimensions :cite:`GoerzEPJQT2015,GoerzNPJQI17,CuiQST17,PatschPRA18`
and open quantum systems :cite:`GoerzNJP2014` in particular
require considerable numerical effort to optimize. Compared to the
Fortran and C/C++ languages traditionally used for scientific computing,
and more recently Julia :cite:`BezansonSIREV2017`, pure
Python code usually performs slower by two to three orders of
magnitude :cite:`AkeretAC2015,EichhornCSJ2018`. Thus, for
hard optimization problems that require several thousand iterations to
converge, the Python implementation provided by the :mod:`krotov` package
may not be sufficiently fast. In this case, it may be desirable to
implement the entire optimization and time propagation in a single, more
efficient (compiled) language. Our Python implementation of Krotov's
method puts an emphasis on clarity, and the documentation provides
detailed explanations of all necessary concepts, especially the correct
:ref:`TimeDiscretization` and the possibility to parallelize the optimization.
Thus, the :mod:`krotov` package can serve as a reference implementation,
leveraging Python's reputation as "executable pseudocode", and as a foundation
against which to test other implementations.

.. _QuTiP library: http://qutip.org
.. _Jupyter notebook: https://jupyter.org
