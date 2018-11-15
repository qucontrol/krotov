.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

Report Bugs
-----------

Report bugs at https://github.com/qucontrol/krotov/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.


Submit Feedback
---------------

The best way to send feedback is to file an issue at https://github.com/qucontrol/krotov/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)


Fix Bugs / Implement Features
-----------------------------

Look through the GitHub issues for bugs or feature requests. Anybody is welcome to submit a pull request for open issues.


Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in ``features.rst``.
3. Check https://travis-ci.org/qucontrol/krotov/pull_requests
   and make sure that the tests pass for all supported Python versions.


Get Started!
------------

Ready to contribute? Follow `Aaron Meurer's Git Workflow Notes`_ (with ``qucontrol/krotov`` instead of ``sympy/sympy``)

In short, if you are not a member of the `qucontrol organization`_,

1. Clone the repository from ``git@github.com:qucontrol/krotov.git``
2. Fork the repo on GitHub to your personal account.
3. Add your fork as a remote.
4. Pull in the latest changes from the master branch.
5. Create a topic branch
6. Make your changes and commit them (testing locally)
7. Push changes to the topic branch on *your* remote
8. Make a pull request against the base master branch through the Github website of your fork.

The project contains a ``Makefile`` to help with development tasts. In your checked-out clone, do

.. code-block:: console

    $ make help

to see the available make targets.

If you are a member of the `qucontrol organization`_, there is no need to fork
``krotov``: you can directly pull and push to ``git@github.com:qucontrol/krotov.git``.

It is strongly recommended that you use the conda_ package manager. The
``Makefile`` relies on conda to create local testing and documentation building
environements (``make test`` and ``make docs``).

Alternatively, you may  use ``make develop-test`` and ``make develop-docs`` to
run the tests or generate the documentation within your active Python
environment. You will have to ensure that all the necessary dependencies are
installed. Also, you will not be able to test the package against all supported
Python versions.
You still can (and should) look at https://travis-ci.org/qucontrol/krotov/ to check that your commits pass all tests.


.. _conda: https://conda.io/docs/

.. _Aaron Meurer's Git Workflow Notes:  https://www.asmeurer.com/git-workflow/

.. _qucontrol organization: https://github.com/qucontrol


Testing
-------

The Krotov package includes a full test-suite using pytest_. We strive for a `test coverage`_ above 90%.

From a checkout of the ``krotov`` repository, assuming conda_ is installed, you can use

.. code-block:: console

    $ make test

to run the entire test suite.

The tests are organized in the ``tests`` subfolder. It includes python scripts
whose name start with ``test_``, which contain functions whose names also start
with ``test_``. Any such functions in any such files are picked up by `pytest`_
for testing. In addition, doctests_ from any docstring or any documentation
file (``*.rst``) are picked up (by the `pytest doctest plugin`_). Lastly, all
`example notebooks <Contribute Examples>`_ are validated as a test, through
the `nbval plugin`_.

.. _test coverage: https://coveralls.io/github/qucontrol/krotov?branch=master
.. _pytest: https://docs.pytest.org/en/latest/
.. _doctests: https://docs.python.org/3.7/library/doctest.html
.. _pytest doctest plugin: https://docs.pytest.org/en/latest/doctest.html
.. _nbval plugin: https://nbval.readthedocs.io/en/latest/


Write Documentation
-------------------

The ``krotov`` package could always use more documentation, whether
as part of the official docs, in docstrings, or even on the web in blog posts,
articles, and such.

The package documentation is generated with Sphinx_, the
documentation (and docstrings) are formatted using the
`Restructured Text markup language`_ (file extension ``rst``).
See also the `Matplotlib Sphinx Sheet sheet`_ for some helpful tips.

Each function or class must have a docstring_; this docstring must
be written in the `"Google Style" format`_ (as implemented by
Sphinx' `napoleon extension`_). Docstrings and any other part of the
documentation can include `mathematical formulas in LaTeX syntax`_
(using mathjax_). In addition to Sphinx' normal syntax for inline math
(``:math:`x```), you may also use easier-to-read dollar signs (``$x$``).
The Krotov package defines some custom tex macros for quantum mechanics, which
you are strongly encouraged to use. These include:

* ``\bra``, e.g. ``$\bra{\Psi}$`` for :math:`\bra{\Psi}` (or ``\\Bra{}`` for auto-resizing).
  Do not use ``\langle``/``\rangle``/``\vert`` manually!
* ``\ket``, e.g. ``$\ket{\Psi}$`` for :math:`\ket{\Psi}` (or ``\Ket{}`` for auto-resizing)
* ``\Braket``, e.g. ``$\Braket{\Phi}{\Psi}$`` for :math:`\Braket{\Phi}{\Psi}`.
* ``\Op`` for for quantum operators, e.g. ``$\Op{H}$`` for :math:`\Op{H}`.
* ``\Abs`` for absolute values, e.g. ``$\Abs{x}$`` for :math:`\Abs{x}`.
* ``\AbsSq``  for the absolute-square, e.g. ``$\AbsSq{\Braket{\Phi}{\Psi}}$`` for :math:`\AbsSq{\Braket{\Phi}{\Psi}}`
* ``\Norm`` for the norm, e.g. ``$\Norm{\ket{\Psi}}$`` for :math:`\Norm{\ket{\Psi}}`
* ``\identity`` for the identity operator, :math:`\identity`
* ``\Liouville`` for the Liouvillian symbol, :math:`\Liouville`

You may use the BibTeX_ plugin for citations.

At any point, from a checkout of the ``krotov`` repository (and
assuming you have conda_ installed), you may run

.. code-block:: console

    $ make docs

to generate the documentation locally.

.. _Sphinx: http://www.sphinx-doc.org/en/master/
.. _Restructured Text markup language: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
.. _docstring: https://www.python.org/dev/peps/pep-0257/
.. _"Google Style" format: http://www.sphinx-doc.org/en/master/usage/extensions/example_google.html#example-google
.. _napoleon extension: http://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
.. _mathematical formulas in LaTeX syntax: http://www.sphinx-doc.org/en/1.6/ext/math.html
.. _mathjax: http://www.sphinx-doc.org/en/master/usage/extensions/math.html#module-sphinx.ext.mathjax
.. _BibTeX: https://sphinxcontrib-bibtex.readthedocs.io/en/latest/
.. _Matplotlib Sphinx Sheet sheet: https://matplotlib.org/sampledoc/cheatsheet.html


Contribute Examples
-------------------

Examples should be contributed in the form of `Jupyter notebooks`_.

.. _Jupyter notebooks: https://jupyter.readthedocs.io/en/latest/index.html

Example notebooks are automatically rendered as part of the documentation
(:ref:`krotov-example-notebooks`), and they are also verified by the automated
tests. For this to work properly, the following steps must be taken:

* Put all imports near the top of the notebook, with ``# NBVAL_IGNORE_OUTPUT``
  as the first line. Use the `watermark`_ package to print out the versions of
  imported packages. For example::

    # NBVAL_IGNORE_OUTPUT
    %load_ext watermark
    import qutip
    import numpy as np
    import scipy
    import matplotlib
    import matplotlib.pylab as plt
    %watermark -v --iversions

* Put the notebook in the folder ``docs/notebooks/``.

* Before committing, re-evaluate all example notebooks in a well-defined
  virtual environment by running

    .. code-block:: console

        $ make notebooks

* Check that the examples can be verified across different Python version by running

    .. code-block:: console

        $ make test

* You may also verify that the example is properly integrated in the documentation by running

    .. code-block:: console

        $ make docs

.. _watermark: https://github.com/rasbt/watermark



Developers' How-tos
-------------------

The following assumes your current working directory is a checkout of
``krotov``, and that you have successfully run ``make test`` (which creates
some local virtual environments that development relies on).

* **How to reference a Github issue in a commit message**

    Simply put e.g. ``#14`` anywhere in your commit message, and Github will
    automatically link to your commit on the page for issue number 14.

    You may also use something like ``Closes #14`` as the last line of your
    commit message to automatically close the issue.
    See `Closing issues using keywords`_ for details.

    Also note the general `Commit Message Guidelines`_.


* **How to run a subset of tests**

    To run e.g. only the tests defined in ``tests/test_krotov.py``, use::

        $ ./.venv/py36/bin/pytest tests/test_krotov.py

    See the `pytest test selection docs`_ for details.

* **How to run only as single test**

    Decorate the test with e.g. ``@pytest.mark.xxx``, and then run, e.g::

        $ ./.venv/py36/bin/pytest -m xxx tests/

    See the `pytest documentation on markers`_ for details.

* **How to run only the doctests**

    Run the following::

    $ ./.venv/py36/bin/pytest --doctest-modules src

* **How to go into an interactive debugger.**

    Optionally, install the `pdbpp` package into the testing environment, for a
    better experience::

        $ ./.venv/py36/bin/python -m pip install pdbpp

    Then:

    - before the line where you went to enter the debugger, insert a line::

        from IPython.terminal.debugger import set_trace; set_trace() # DEBUG

    - Run ``pytest`` with the option ``-s``, e.g.::

        $ ./.venv/py36/bin/pytest -m xxx -s tests/

    You may also see the `pytest documentation on automatic debugging`_.


.. _Closing issues using keywords: https://help.github.com/articles/closing-issues-using-keywords/
.. _Commit Message Guidelines: https://gist.github.com/robertpainsi/b632364184e70900af4ab688decf6f53
.. _pytest test selection docs: https://docs.pytest.org/en/latest/usage.html#specifying-tests-selecting-tests
.. _pytest documentation on markers: https://docs.pytest.org/en/latest/example/markers.html
.. _pytest documentation on automatic debugging: https://docs.pytest.org/en/latest/usage.html#dropping-to-pdb-python-debugger-on-failures
