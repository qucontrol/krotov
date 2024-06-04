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

* `Juqbox.jl <https://github.com/LLNL/Juqbox.jl>`_ (Julia) :cite:`AndersPeterssonSJSC2022` -- Quantum optimal control framework implementing methods developed at LLNL (symplectic time integration, parameterization of the control functions using B-splines with carrier waves)

* `Quandary <https://github.com/LLNL/quandary>`_ (C++) -- Optimal control for open and closed quantum systems for controls as B-spline polynomials

* `C3 <https://github.com/q-optimize/c3>`_ (Python) :cite:`WittlerPRA2021`  -- An integrated tool-set for Control, Calibration and Characterization

* `QuantumControl.jl <https://github.com/JuliaQuantumControl/QuantumControl.jl>`_ (Julia) -- Comprehensive Framework for Quantum Control, within the the `JuliaQuantumControl org <https://github.com/JuliaQuantumControl>`_

* `Krotov.jl <https://github.com/JuliaQuantumControl/Krotov.jl>`_ (Julia) -- Reimplementation of this package in Julia, with further features such as automatic differentiation. Part of the `JuliaQuantumControl org <https://github.com/JuliaQuantumControl>`_.

* `QControl.jl <https://github.com/Phionx/QControl.jl>`_ (Julia) -- Quantum Control via Constrained Trajectory Optimization

* `qontrol <https://github.com/dkweiss31/qontrol>`_ (Python) -- Quantum optimal control package built on top of `dynamics <https://github.com/dynamiqs/dynamiqs>`_, using JAX

* `PRONTO.jl <https://github.com/narijauskas/PRONTO.jl>`_ (Julia) --  PRojection-Operator-Based Newton’s Method for Trajectory Optimization

* `Piccolo.jl <https://github.com/kestrelquantum/Piccolo.jl>`_ (Julia) -- Quantum optimal control using the Pade Integrator Collocation (Piccolo) method


Also:

* `Pulser <https://github.com/pasqal-io/Pulser>`_ (Python) :cite:`SilverioQ2022` -- Framework for composing, simulating and executing pulse sequences for neutral-atom quantum devices

* `Bloqade.jl <https://github.com/QuEraComputing/Bloqade.jl>`_ (Julia) -- Quantum dynamics on neutral-atom architectures

* `Mitiq <https://github.com/unitaryfund/mitiq>`_ (Python) -- Toolkit for implementing error mitigation techniques on quantum computers

* `QuantumOptics.jl <https://github.com/qojulia/QuantumOptics.jl>`_ (Julia) -- Framework to simulate various kinds of open quantum systems, inspired by QuTiP

* `QuantumToolbox.jl <https://github.com/qutip/QuantumToolbox.jl>`_ (Julia) -- Translation of QuTiP to Julia, aiming to preserve the QuTiP API as closely as possible

* `QuantumSavory.jl <http://qs.quantumsavory.org/stable/>`_ (Julia) -- Multi-formalism simulator for noisy quantum communication and computation hardware

* `QuanEstimation <https://github.com/QuanEstimation/QuanEstimation>`_ (Python/Julia) -- Toolkit for quantum parameter estimation


Accessories
-----------

The following packages integrate closely with :mod:`krotov`.

* `weylchamber <https://github.com/qucontrol/weylchamber>`_ (Python) -- Package for analyzing two-qubit gates in the Weyl chamber. Provides :ref:`Local-Invariants <HowtoLIOptimization>` and :ref:`Perfect Entangler <HowtoPEOptimization>` functionals for use with :mod:`krotov`.
