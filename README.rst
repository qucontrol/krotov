=====================
Krotov Python Package
=====================
.. image:: https://img.shields.io/pypi/v/krotov.svg
   :alt: Krotov on the Python Package Index
   :target: https://pypi.python.org/pypi/krotov
.. image:: https://badges.gitter.im/qucontrol_krotov/Lobby.svg
   :alt: Join the chat at https://gitter.im/qucontrol_krotov/Lobby
   :target: https://gitter.im/qucontrol_krotov/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
.. image:: https://img.shields.io/travis/qucontrol/krotov.svg
   :alt: Travis Continuous Integration
   :target: https://travis-ci.org/qucontrol/krotov
.. image:: https://coveralls.io/repos/github/qucontrol/krotov/badge.svg?branch=master&service=github:alt: Coveralls
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

This implementation is built on top of `QuTiP`_.

Development of the ``krotov`` package happens on `Github`_. You can read the full documentation at `ReadTheDocs`_.


.. _QDYN Fortran library: https://www.qdyn-library.net
.. _QuTiP: http://qutip.org
.. _ReadTheDocs: https://krotov.readthedocs.io/en/latest/


Installation
------------
To install the latest released version of ``krotov``, run this command in your terminal:

.. code-block:: console

    $ pip install krotov

This is the preferred method to install the ``krotov`` package, as it will always install the most recent stable release.

**Note: the latest released version is a placeholder release that is non-functional**

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


To install the latest development version of ``krotov`` from `Github`_.

.. code-block:: console

    $ pip install git+https://github.com/qucontrol/krotov.git@master#egg=krotov

.. _Github: https://github.com/qucontrol/krotov

Usage
-----

To use the ``krotov`` package in a Python project::

    import krotov


