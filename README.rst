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
.. image:: https://img.shields.io/badge/License-BSD-green.svg
   :alt: BSD License
   :target: https://opensource.org/licenses/BSD-3-Clause
.. image:: https://readthedocs.org/projects/krotov/badge/?version=latest
   :alt: Documentation Status
   :target: https://krotov.readthedocs.io/en/latest/?badge=latest
.. image:: https://img.shields.io/badge/arXiv-1902.11284-red.svg
   :alt: arXiv
   :target: https://arxiv.org/abs/1902.11284

Python implementation of Krotov's method for quantum optimal control.

This implementation follows the original implementation in the `QDYN Fortran library`_.
The method is described in detail in `D. M. Reich, M. Ndong, and C. P. Koch, J. Chem. Phys. 136, 104103 (2012) <https://doi.org/10.1063/1.3691827>`_ (`arXiv:1008.5126 <http://arxiv.org/abs/1008.5126>`_).

The ``krotov`` package is built on top of `QuTiP`_.

Development happens on `Github`_. You can read the full documentation at `ReadTheDocs`_.

If you use the ``krotov`` package in your research, please `cite it <https://krotov.readthedocs.io/en/stable/01_overview.html#citing-the-krotov-package>`_.

.. _QDYN Fortran library: https://www.qdyn-library.net
.. _QuTiP: http://qutip.org
.. _ReadTheDocs: https://krotov.readthedocs.io/en/stable/


Purpose
-------

Optimal control is a cornerstones of quantum technology: relying not
just on a passive understanding of quantum mechanics, but on the *active*
utilization of the quantum properties of matter. Quantum optimal control asks
how to manipulate the dynamics of a quantum system to behave in some desired
way. This is essential for the realization of quantum computers, and
related technologies such as quantum sensing.  See e.g. `Glaser et al. Eur.
Phys. J. D 69, 279 (2015)`_ for an overview of methods, applications, and
current research directions. Quantum technology and thus quantum control are
the focus of several large-scale national and super-national research
endeavors, such as the U.S. $1 billion `National Quantum Initiative Act`_ and
the €1 billion `European Quantum Flagship program`_.

Krotov's method is one of the two leading gradient-based optimization
algorithms used in numerical quantum optimal control. It simulates the dynamics
of a quantum system under a set of initial controls, and evaluates the
result with respect to an optimization functional to be minimized. It then
iteratively modifies the controls to guarantee a monotonically decreasing value
in the optimization functional. To date, there has not been an open source
implementation of the method. This package provides that missing
implementation.

The choice of Python as an implementation language is due to Python's easy to learn
syntax, expressiveness, and immense popularity in the scientific community.
Moreover, the `QuTiP`_ library exists to provide the foundations of
numerically describing quantum systems, and already includes basic versions of
some of the other popular algorithms in quantum control, the gradient-based
GRAPE and the gradient-free CRAB. The availability of the `Jupyter notebook`_
framework provides an ideal platform for showing the use of the method.

The Krotov package targets both students wishing to enter
the field of quantum control, and researchers in the field. By providing a
comprehensive set of examples, we enable users of our package to
explore the formulation of typical control problems, and to understand how
Krotov's method can solve them. These examples are inspired by
recent publications, and thus show the use of the method at the cutting edge of
research. Optimal control is also increasingly important in the design of
experiments, and we hope that the availability of an easy to use implementation
of Krotov's will facilitate this further.

The use of Python implies that for large-scale control problems, performance
may become a significant issue. In this case, it may be necessary to implement
Krotov's method in a more efficient (compiled) language. While the method as
such is relatively straightforward, there are some subtleties involved. Our
implementation puts an emphasis on clarity, and the documentation provides
detailed explanations of all necessary concepts.  Thus, the Krotov package can
serve as a reference implementation, leveraging Python's reputation as
"executable pseudocode", and as a foundation against which to test other
implementations.

.. _Glaser et al. Eur. Phys. J. D 69, 279 (2015): https://link.springer.com/article/10.1140%2Fepjd%2Fe2015-60464-1
.. _European Quantum Flagship program: https://qt.eu/about/
.. _National Quantum Initiative Act: https://www.forbes.com/sites/alexknapp/2018/12/20/congress-just-passed-a-bill-to-accelerate-quantum-computing-heres-what-it-does/#20b5d2c22ef8


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

If you don't have `pip`_ installed, the `Python installation guide`_, respectively the `Python Packaging User Guide`_ can guide
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
  |krotov.Objective|_
* call |krotov.optimize_pulses|_ to perform an optimization of an arbitrary
  number of control fields over all the objectives.

.. |krotov.Objective| replace:: ``krotov.Objective``
.. _krotov.Objective: https://krotov.readthedocs.io/en/stable/API/krotov.objectives.html#krotov.objectives.Objective

.. |krotov.optimize_pulses| replace:: ``krotov.optimize_pulses``
.. _krotov.optimize_pulses: https://krotov.readthedocs.io/en/stable/API/krotov.optimize.html#krotov.optimize.optimize_pulses

See `Using Krotov with QuTiP <https://krotov.readthedocs.io/en/stable/07_qutip_usage.html#using-krotov-with-qutip>`_ and `Examples <https://krotov.readthedocs.io/en/stable/08_examples.html>`_ for details.

.. _Jupyter notebook: http://jupyter.org
