# Krotov Python Package

[![Source code on Github](https://img.shields.io/badge/github-qucontrol/krotov-blue.svg)](https://github.com/qucontrol/krotov)
[![Documentation](https://img.shields.io/badge/docs-gh--pages-blue.svg)](https://qucontrol.github.io/krotov)
[![Krotov on the Python Package Index](https://img.shields.io/pypi/v/krotov.svg)](https://pypi.python.org/pypi/krotov)
[![Docs](https://github.com/qucontrol/krotov/actions/workflows/docs.yml/badge.svg?branch=master)](https://github.com/qucontrol/krotov/actions?query=workflow%3ADocs)
[![Tests](https://github.com/qucontrol/krotov/actions/workflows/test.yml/badge.svg?branch=master)](https://github.com/qucontrol/krotov/actions?query=workflow%3ATests)
[![Codecov](https://codecov.io/gh/qucontrol/krotov/branch/master/graph/badge.svg)](https://codecov.io/gh/qucontrol/krotov)
[![BSD License](https://img.shields.io/badge/License-BSD-green.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/qucontrol/krotov/v1.3.0?filepath=docs%2Fnotebooks)
[![DOI](https://img.shields.io/badge/DOI-10.21468/SciPostPhys.7.6.080-blue.svg)](https://doi.org/10.21468/SciPostPhys.7.6.080)

Python implementation of Krotov's method for quantum optimal control.

This implementation follows the original implementation in the [QDYN
Fortran library](https://www.qdyn-library.net).

The `krotov` package is built on top of [QuTiP](http://qutip.org).

Development happens on [Github](https://github.com/qucontrol/krotov).
You can read the full documentation
[online](https://qucontrol.github.io/krotov).

If you use the `krotov` package in your research, please [cite
it](https://qucontrol.github.io/krotov/v1.3.0/01_overview.html#citing-the-krotov-package).

## Purpose

Optimal control is a cornerstone of quantum technology: relying not just
on a passive understanding of quantum mechanics, but on the *active*
utilization of the quantum properties of matter. Quantum optimal control
asks how to manipulate the dynamics of a quantum system in some desired
way. This is essential for the realization of quantum computers and
related technologies such as quantum sensing.

Krotov's method and GRAPE are the two leading gradient-based
optimization algorithms used in numerical quantum optimal control.
Krotov's method distinguishes itself by guaranteeing monotonic
convergence for near-continuous control fields. This makes is
particularly useful for exploring the limits of controllability in a
physical system. While GRAPE is found in various software packages,
there has not been an open source implementation of Krotov's method to
date. Our package provides that missing implementation.

The Krotov package targets both students wishing to enter the field of
quantum control and researchers in the field. It was designed towards
the following goals:

- Leverage the [QuTiP](http://qutip.org) library as a platform for
  numerically describing quantum systems.
- Provide a collection of examples inspired by recent publications in
  the [Jupyter notebook](https://jupyter.org) format, allowing for
  interactive exploration of the method.
- Define a general interface for formulating *any* quantum control
  problem, which may extend to other optimization methods in the future.
- Serve as a reference implementation of Krotov's method, and as a
  foundation against which to test other implementations.
- Enable the more widespread use of Krotov's method, for example in the
  design of experiments.

## Further Information

For further information, including installation and usage instructions, see the
documentation at https://qucontrol.github.io/krotov.
