=====================
Krotov Python Package
=====================
.. image:: https://img.shields.io/badge/github-qucontrol/krotov-blue.svg
   :alt: Source code on Github
   :target: https://github.com/qucontrol/krotov
.. image:: https://img.shields.io/pypi/v/krotov.svg
   :alt: Krotov on the Python Package Index
   :target: https://pypi.python.org/pypi/krotov
.. image:: https://badges.gitter.im/qucontrol_krotov/Lobby.svg
   :alt: Join the chat at https://gitter.im/qucontrol_krotov/Lobby
   :target: https://gitter.im/qucontrol_krotov/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
.. image:: https://img.shields.io/travis/qucontrol/krotov.svg
   :alt: Travis Continuous Integration
   :target: https://travis-ci.org/qucontrol/krotov
.. image:: https://ci.appveyor.com/api/projects/status/1cbm24w04jmxjpjh?svg=true
   :alt: AppVeyor Continuous Integration
   :target: https://ci.appveyor.com/project/goerz/krotov
.. image:: https://img.shields.io/coveralls/github/qucontrol/krotov/master.svg
   :alt: Coveralls
   :target: https://coveralls.io/github/qucontrol/krotov?branch=master
.. image:: https://readthedocs.org/projects/krotov/badge/?version=latest
   :alt: Documentation Status
   :target: https://krotov.readthedocs.io/en/latest/?badge=latest
.. image:: https://img.shields.io/badge/License-BSD-green.svg
   :alt: BSD License
   :target: https://opensource.org/licenses/BSD-3-Clause

Python implementation of Krotov's method for quantum optimal control.

This implementation follows the original implementation in the `QDYN Fortran library`_.
The method is described in detail in `D. M. Reich, M. Ndong, and C. P. Koch, J. Chem. Phys. 136, 104103 (2012) <https://doi.org/10.1063/1.3691827>`_ (`arXiv:1008.5126 <http://arxiv.org/abs/1008.5126>`_)

The ``krotov`` package is built on top of `QuTiP`_.

Development happens on `Github`_. You can read the full documentation at `ReadTheDocs`_.


.. _QDYN Fortran library: https://www.qdyn-library.net
.. _QuTiP: http://qutip.org
.. _ReadTheDocs: https://krotov.readthedocs.io/en/latest/


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

.. code-block:: console

    $ conda create -n qucontrolenv python=3.6
    $ conda activate qucontrolenv
    $ conda config --append channels conda-forge
    $ conda install qutip

.. _core packages of the Python scientific ecosystem: https://www.scipy.org
.. _QuTiP's installation instructions: http://qutip.org/docs/latest/installation.html
.. _virtual environment: https://docs.python.org/3/glossary.html#term-virtual-environment
.. _wheels: https://packaging.python.org/tutorials/installing-packages/#source-distributions-vs-wheels
.. _QuTiP currently does not provide wheels: https://github.com/qutip/qutip/issues/933
.. _conda: https://conda.io/docs/index.html
.. _Miniconda: https://conda.io/miniconda.html


Installation
------------
To install the latest released version of ``krotov`` into your current (conda)
environment, run this command in your terminal:

.. code-block:: console

    $ pip install krotov

This is the preferred method to install the ``krotov`` package, as it will always install the most recent stable release.

You may also do

.. code-block:: console

    $ pip install krotov[dev,extras]

to install additional development dependencies, including packages required to run the example notebooks.

If you don't have `pip`_ installed, this `Python installation guide`_, respectively the `Python Packaging User Guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _Python Packaging User Guide: https://packaging.python.org/tutorials/installing-packages/


To install the latest development version of ``krotov`` from `Github`_:

.. code-block:: console

    $ pip install git+https://github.com/qucontrol/krotov.git@master#egg=krotov

.. _Github: https://github.com/qucontrol/krotov

Usage
-----

To use Krotov's method for quantum optimal control in a Python script or
`Jupyter notebook`_, start with::

    import krotov

Then,

* define the necessary quantum operators and states using `QuTiP`_.
* create a list of objectives, as instances of
  :class:`krotov.Objective <krotov.objectives.Objective>`
* call :func:`krotov.optimize_pulses <krotov.optimize.optimize_pulses>` to
  perform an optimization of an arbitrary number of control fields over all the
  objectives.

See :ref:`using-krotov-with-qutip` and :ref:`krotov-example-notebooks` for details.

.. _Jupyter notebook: http://jupyter.org
