Related Software
================

Is there any software missing from this list? `Open an issue! <https://github.com/qucontrol/krotov/issues/new?assignees=&labels=documentation&template=related-software.md&title=>`_


Other implementations of quantum control
----------------------------------------

* `QDYN <https://www.qdyn-library.net>`_ (Fortran) -- comprehensive package for quantum dynamics and control. The implementation of the :mod:`krotov` package was inspired by QDYN.

* `QuTiP <http://qutip.org>`_ (Python) -- Quantum Toolbox in Python :cite:`JohanssonCPC2012,JohanssonCPC2013`. Provides the basic data structures for the :mod:`krotov`. Contains implementation of :ref:`GRAPE <GrapeInQutip>` and :ref:`CRAB <GradientFreeOptimization>`.

* `WavePacket <https://sourceforge.net/p/wavepacket/wiki/Home/>`_ (Matlab) :cite:`SchmidtCPC2018` -- Implemenentation of :ref:`monotonically converging iterative schemes <IterativeSchemes>` for wavepacket dynamics

* `DYNAMO <https://github.com/shaimach/Dynamo>`_ (Matlab) :cite:`MachnesPRA11` -- Implementation of :ref:`GRAPE`, including concurrent/sequenetial/hybrid schemes

* `QEngine <https://gitlab.com/quatomic/qengine>`_ (C++) :cite:`SorensenCPC2019` -- Implementation of :ref:`"GRadient Optimization Using Parametrization" (GROUP) <GRAPE>` :cite:`SorensenPRA2018`

Accessories
-----------

The following packages integrate closely with :mod:`krotov`.

* `weylchamber <https://github.com/qucontrol/weylchamber>`_ (Python) -- Package for analyzing two-qubit gates in the Weyl chamber. Provides :ref:`Local-Invariants <HowtoLIOptimization>` and :ref:`Perfect Entangler <HowtoPEOptimization>` functionals for use with :mod:`krotov`.
