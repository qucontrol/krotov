=====================
Krotov Python Package
=====================

.. image:: https://img.shields.io/badge/github-qucontrol/krotov-blue.svg
   :alt: Source code on Github
   :target: https://github.com/qucontrol/krotov
.. image:: https://img.shields.io/badge/docs-gh--pages-blue.svg
   :alt: Documentation
   :target: https://qucontrol.github.io/krotov
.. image:: https://img.shields.io/pypi/v/krotov.svg
   :alt: Krotov on the Python Package Index
   :target: https://pypi.python.org/pypi/krotov
.. image:: https://badges.gitter.im/qucontrol_krotov/Lobby.svg
   :alt: Join the chat at https://gitter.im/qucontrol_krotov/Lobby
   :target: https://gitter.im/qucontrol_krotov/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
.. image:: https://github.com/qucontrol/krotov/workflows/Docs/badge.svg?branch=master
   :alt: Docs
   :target: https://github.com/qucontrol/krotov/actions?query=workflow%3ADocs
.. image:: https://github.com/qucontrol/krotov/workflows/Tests/badge.svg?branch=master
   :alt: Tests
   :target: https://github.com/qucontrol/krotov/actions?query=workflow%3ATests
.. image:: https://codecov.io/gh/qucontrol/krotov/branch/master/graph/badge.svg
   :alt: Codecov
   :target: https://codecov.io/gh/qucontrol/krotov
.. image:: https://img.shields.io/badge/License-BSD-green.svg
   :alt: BSD License
   :target: https://opensource.org/licenses/BSD-3-Clause
.. image:: https://mybinder.org/badge_logo.svg
   :alt: Launch Binder
   :target: https://mybinder.org/v2/gh/qucontrol/krotov/v1.2.1?filepath=docs%2Fnotebooks
.. image:: https://img.shields.io/badge/DOI-10.21468/SciPostPhys.7.6.080-blue.svg
   :alt: DOI
   :target: https://doi.org/10.21468/SciPostPhys.7.6.080

Python implementation of Krotov's method for quantum optimal control.

This implementation follows the original implementation in the `QDYN Fortran library`_.

The :mod:`krotov` package is built on top of `QuTiP`_.

Development happens on `Github`_. You can read the full documentation `online`__ or `download a PDF version`_.

.. _Documentation: https://qucontrol.github.io/krotov
__ Documentation_

If you use the :mod:`krotov` package in your research, please :ref:`cite it <CitingKrotov>`.

.. _QDYN Fortran library: https://www.qdyn-library.net
.. _QuTiP: http://qutip.org
.. _download a PDF version: https://github.com/qucontrol/krotov/tree/master/docs/pdf


Purpose
-------

Optimal control is a cornerstone of quantum technology: relying not
just on a passive understanding of quantum mechanics, but on the *active*
utilization of the quantum properties of matter. Quantum optimal control asks
how to manipulate the dynamics of a quantum system in some desired
way. This is essential for the realization of quantum computers and
related technologies such as quantum sensing.

Krotov's method and GRAPE are the two leading gradient-based optimization
algorithms used in numerical quantum optimal control. Krotov's method
distinguishes itself by guaranteeing monotonic convergence for near-continuous
control fields. This makes is particularly useful for exploring the limits of
controllability in a physical system.
While GRAPE is found in various software packages, there has not been an open
source implementation of Krotov's method to date. Our package provides that
missing implementation.

The Krotov package targets both students wishing to enter the field
of quantum control and researchers in the field. It was designed towards
the following goals:

* Leverage the `QuTiP`_ library as a platform for numerically describing
  quantum systems.
* Provide a collection of examples inspired by recent publications in
  the `Jupyter notebook`_ format, allowing for interactive exploration of the
  method.
* Define a general interface for formulating *any* quantum control problem,
  which may extend to other optimization methods in the future.
* Serve as a reference implementation of Krotov's method, and as a foundation
  against which to test other implementations.
* Enable the more widespread use of Krotov's method, for example in the design
  of experiments.


Prerequisites
-------------

The Krotov package is available for Python versions >= 3.5. Its main dependency is `QuTiP`_
(apart from the `core packages of the Python scientific ecosystem`_).
Thus, you should consider `QuTiP's installation instructions`_.

In any case, using some sort of `virtual environment`_ is strongly encouraged.
Most packages in the Python scientific ecosystem are now available as
`wheels`_, making installation via `pip`_ easy. However, `QuTiP currently does
not provide wheels`_. Thus, on systems that do not have the necessary compilers
installed (Windows, macOS), the `conda`_ package manager provides a good solution.

Assuming ``conda`` is installed (e.g. through `Miniconda`_), the following
commands set up a virtual (conda) environment into which the Krotov package can
then be installed:

.. code-block:: shell

    conda create -n qucontrolenv "python=3.7"
    conda activate qucontrolenv
    conda config --append channels conda-forge
    conda install qutip

.. _core packages of the Python scientific ecosystem: https://www.scipy.org
.. _QuTiP's installation instructions: http://qutip.org/docs/latest/installation.html
.. _virtual environment: https://docs.python.org/3/glossary.html#term-virtual-environment
.. _wheels: https://packaging.python.org/tutorials/installing-packages/#source-distributions-vs-wheels
.. _QuTiP currently does not provide wheels: https://github.com/qutip/qutip/issues/933
.. _conda: https://conda.io/docs/index.html
.. _Miniconda: https://conda.io/miniconda.html


Installation
------------
To install the latest released version of :mod:`krotov` into your current (conda)
environment, run this command in your terminal:

.. code-block:: shell

    python -m pip install krotov

This is the preferred method to install the :mod:`krotov` package, as it will always install the most recent stable release.

You may also do

.. code-block:: shell

    python -m pip install krotov[dev,extras]

to install additional development dependencies, including packages required to run the example notebooks.

If you don't have `pip`_ installed, the `Python installation guide`_, respectively the `Python Packaging User Guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _Python Packaging User Guide: https://packaging.python.org/tutorials/installing-packages/


To install the latest development version of :mod:`krotov` from `Github`_:

.. code-block:: shell

    python -m pip install git+https://github.com/qucontrol/krotov.git@master#egg=krotov

.. _Github: https://github.com/qucontrol/krotov

Usage
-----

To use Krotov's method for quantum optimal control in a Python script or
`Jupyter notebook`_, start with::

    import krotov
    import qutip

Then,

1. define the necessary quantum operators and states using `QuTiP`_.
2. create a list of objectives, as instances of :class:`krotov.Objective <krotov.objectives.Objective>`.
3. call :func:`krotov.optimize_pulses <krotov.optimize.optimize_pulses>` to perform an optimization of an arbitrary
   number of control fields over all the objectives.


See :ref:`using-krotov-with-qutip` and :ref:`krotov-example-notebooks` for details.

.. _Jupyter notebook: https://jupyter.org
