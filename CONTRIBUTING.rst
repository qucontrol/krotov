.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.


Code of Conduct
---------------

Everyone interacting in the ``krotov`` project's code base,
issue tracker, and any communication channels is expected to follow the
`PyPA Code of Conduct`_.

.. _`PyPA Code of Conduct`: https://www.pypa.io/en/latest/code-of-conduct/


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


Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in ``docs/04_features.rst`` and/or ``HISTORY.rst``.
3. Check https://github.com/qucontrol/krotov/actions
   and make sure that the tests pass for all supported Python versions.


Get Started!
------------

Ready to contribute? Follow `Aaron Meurer's Git Workflow Notes`_ (with ``qucontrol/krotov`` instead of ``sympy/sympy``)

In short, if you are not a member of the `qucontrol organization`_,

1. Clone the repository from ``git@github.com:qucontrol/krotov.git``
2. Fork the repo on GitHub to your personal account.
3. Add your fork as a remote.
4. Pull in the latest changes from the master branch.
5. Create a topic branch.
6. Make your changes and commit them (testing locally).
7. Push changes to the topic branch on *your* remote.
8. Make a pull request against the base master branch through the Github website of your fork.

The project uses tox_ for automated testing accross multiple versions of Python
and for various development tasks such as linting and generating the
documentation. See :ref:`DevelopmentPrerequisites` for details.

There is also a ``Makefile`` that wraps around tox, for
convenience on Unix-based systems. In your checked-out clone, run

.. code-block:: shell

    make help

to see the available make targets. If you cannot use ``make``, but want to use
``tox`` directly (e.g., on Windows), run

.. code-block:: shell

    tox -av

to see a list of tox environments and a description. For the initial
configuration of tox environments, you may have to run

.. code-block:: shell

    tox -e bootstrap

in order to set up the ``tox.ini`` configuration file.


If you are a member of the `qucontrol organization`_, there is no need to fork
``krotov`` - you can directly pull and push to ``git@github.com:qucontrol/krotov.git``.

.. _tox: https://tox.readthedocs.io

.. _Aaron Meurer's Git Workflow Notes:  https://www.asmeurer.com/git-workflow/

.. _qucontrol organization: https://github.com/qucontrol


.. _DevelopmentPrerequisites:


Development Prerequisites
-------------------------

Contributing to the package's developments requires that you have Python 3.7
and tox_ installed. It is strongly recommended that you also have installations
of all other supported Python versions. The recommended way to install multiple
versions of Python at the same time is through pyenv_ (or pyenv-win_ on
Windows).

Alternatively, you may install conda_ (via the Anaconda_ or Miniconda_
distributions, or also through pyenv_). As ``conda`` can create environments
with any version of Python (independent of which Python version ``conda`` was
originally installed with), this alleviates the need for managing multiple
versions.
The advantage of using conda_ is that you may be able to avoid installing the
compilers necessary for Python extension packages. The disadvantage is that
environment creation is slower and the resulting environments are bigger, and
that you may run into occasional binary incompatibilities between conda packages.

.. warning::
   If you want to use `conda`, you must use the ``tox-conda.ini`` configuration
   file. That is, run all ``make`` comands as e.g.
   ``make TOXINI=tox-conda.ini test`` and ``tox`` commands as e.g.
   ``tox -c tox-conda.ini -e py35-test,py36-test,py37-test``. Alternatively,
   make ``tox-conda.ini`` the default by copying it to ``tox.ini``.

.. _pyenv: https://github.com/pyenv/pyenv
.. _pyenv-win: https://github.com/pyenv-win/pyenv-win
.. _conda: https://conda.io/docs/
.. _Anaconda: https://www.anaconda.com/distribution/
.. _Miniconda: https://conda.io/en/latest/miniconda.html
.. _QuTiP: http://qutip.org


.. _BranchingModel:

Branching Model
---------------

For developers with direct access to the repository,
``krotov`` uses a simple branching model where all
developments happens directly on the ``master`` branch. Releases are tags on
``master``. All commits on ``master`` *should* pass all tests and be
well-documented. This is so that ``git bisect`` can be effective. For any
non-trivial issue, it is recommended to create a topic branch, instead of
working on ``master``. There are no restrictions on commits on topic branches,
they do not need to contain complete documentation, pass any tests, or even be
able to run.

To create a topic-branch named ``issue1``::

    git branch issue1
    git checkout issue1

You can then make commits, and push them to Github to trigger Continuous
Integration testing::

    git push -u origin issue1

Commit early and often! At the same time, try to keep your topic branch
as clean and organized as possible.

* Avoid having a series of meaningless granular commits like "start bugfix",
  "continue development", "add more work on bugfix", "fix typos", and so forth.
  Instead, use ``git commit --amend`` to add to your previous commit. This is
  the ideal way to "commit early and often". You do not have to wait until a
  commit is "perfect"; it is a good idea to make hourly/daily "snapshots" of
  work in progress. Amending a commit also allows you to change the commit
  message of your last commit.
* You can combine multiple existing commits by "squashing" them. For example,
  use ``git rebase -i HEAD~4`` to combined the previous four commits into one.
  See the `"Rewriting History" section of Pro Git book`_ for details (if you
  feel this is too far outside of your git comfort zone, just skip it).
* If you work on a topic branch for a long time, and there is significant work
  on ``master`` in the meantime, periodically rebase your topic branch on the
  current master (``git rebase master``). Avoid merging ``master`` into your
  topic branch. See `Merging vs. Rebasing`_.

If you have already pushed your topic branch to the remote origin, you can
force-push the issue branch (``git push --force``). If you are collaborating
with others on the branch, coordinate with them before force pushing. A
force-push rewrites history. You must never rewrite history on the ``master``
branch (nor will you be able to, as the ``master`` branch is "protected" and
can only be force-pushed to in coordination with the project maintainer).  If
something goes wrong with any advanced "history rewriting", there is always
`"git reflog"`_ as a safety net -- you will never lose work that was committed
before.

When you are done with a topic branch (the issue has been fixed), finish up by
creating a pull-request for merging the branch into ``master`` (follow the
propmts on the Github website).

Summarize the changes of the branch relative to ``master`` in the pull request.


.. _"Rewriting History" section of Pro Git book: https://git-scm.com/book/en/v2/Git-Tools-Rewriting-History
.. _Merging vs. Rebasing: https://www.atlassian.com/git/tutorials/merging-vs-rebasing
.. _"git reflog": https://www.atlassian.com/git/tutorials/rewriting-history/git-reflog


Commit Message Guidelines
-------------------------

Write commit messages according to this template:

.. code-block:: none

    Short (50 chars or less) summary ("subject line")

    More detailed explanatory text. Wrap it to 72 characters. The blank
    line separating the summary from the body is critical (unless you omit
    the body entirely).

    Write your subject line in the imperative: "Fix bug" and not "Fixed
    bug" or "Fixes bug." This convention matches up with commit messages
    generated by commands like git merge and git revert. A properly formed
    git commit subject line should always be able to complete the sentence
    "If applied, this commit will <your subject line here>".

    Further paragraphs come after blank lines.

    - Bullet points are okay, too.
    - Typically a hyphen or asterisk is used for the bullet, followed by a
      single space. Use a hanging indent.

    You should reference any issue that is being addressed in the commit, as
    e.g. "#1" for issue #1. If the commit closes an issue, state this on the
    last line of the message (see below). This will automatically close the
    issue on Github as soon as the commit is pushed there.

    Closes #1

See `Closing issues using keywords`_ for details on references to issues that
Github will understand.


Testing
-------

The Krotov package includes a full test-suite using pytest_. We strive for a `test coverage`_ above 90%.

From a checkout of the ``krotov`` repository, you can use

.. code-block:: shell

    make test

to run the entire test suite, or

.. code-block:: shell

    tox -e py35-test,py36-test,py37-test

if ``make`` is not available.

The tests are organized in the ``tests`` subfolder. It includes python scripts
whose name start with ``test_``, which contain functions whose names also start
with ``test_``. Any such functions in any such files are picked up by `pytest`_
for testing. In addition, doctests_ from any docstring or any documentation
file (``*.rst``) are picked up (by the `pytest doctest plugin`_). Lastly, all
:ref:`example notebooks <ContributeExamples>` are validated as a test, through
the `nbval plugin`_.

.. _test coverage: https://codecov.io/gh/qucontrol/krotov
.. _pytest: https://docs.pytest.org/en/latest/
.. _doctests: https://docs.python.org/3.7/library/doctest.html
.. _pytest doctest plugin: https://docs.pytest.org/en/latest/doctest.html
.. _nbval plugin: https://nbval.readthedocs.io/en/latest/


Code Style
----------

All code must be compatible with :pep:`8`. The line length limit
is 79 characters, although exceptions are permissible if this improves
readability significantly.

Beyond :pep:`8`, this project adopts the `Black code style`_, with
``--skip-string-normalization --line-length 79``. You can
run ``make black-check`` or ``tox -e run-blackcheck`` to check adherence to the
code style, and ``make black`` or ``tox -e run-black`` to apply it.

.. _Black code style: https://github.com/ambv/black/#the-black-code-style

Imports within python modules must be sorted according to the isort_
configuration in ``setup.cfg``. The command ``make isort-check`` or ``tox -e
run-isortcheck`` checks whether all imports are sorted correctly, and ``make
isort`` or ``tox -e run-isort`` modifies all Python modules in-place with the
proper sorting.

.. _isort: https://github.com/timothycrosley/isort#readme

The code style is enforced as part of the test suite, as well as through git
pre-commit hooks that prevent committing code not does not meet the
requirements. These hooks are managed through the `pre-commit framework`_.

.. warning::
   After cloning the ``krotov`` repository, you should run
   ``make bootstrap``, ``tox -e bootstrap``, or ``python scripts/bootstrap.py``
   from within the project root folder. These set up ``tox``, and the
   pre-commit hooks

.. _pre-commit framework: https://pre-commit.com

You may use ``make flake8-check`` or ``tox -e run-flake8`` and ``make
pylint-check`` or ``tox -e run-pylint`` for additional checks on the code with
flake8_ and pylint_, but there is no strict requirement for a perfect score
with either one of these linters. They only serve as a guideline for code that
might be improved.

.. _flake8: http://flake8.pycqa.org
.. _pylint: http://pylint.pycqa.org


.. _write-documentation:

Write Documentation
-------------------

The ``krotov`` package could always use more documentation, whether
as part of the official docs, in docstrings, or even on the web in blog posts,
articles, and such.

The package documentation is generated with Sphinx_, the
documentation (and docstrings) are formatted using the
`Restructured Text markup language`_ (file extension ``rst``).
See also the `Matplotlib Sphinx cheat sheet`_ for some helpful tips.

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
* ``\ket``, e.g. ``$\ket{\Psi}$`` for :math:`\ket{\Psi}` (or ``\Ket{}`` for auto-resizing).
* ``\Braket``, e.g. ``$\Braket{\Phi}{\Psi}$`` for :math:`\Braket{\Phi}{\Psi}`.
* ``\Op`` for quantum operators, e.g. ``$\Op{H}$`` for :math:`\Op{H}`.
* ``\Abs`` for absolute values, e.g. ``$\Abs{x}$`` for :math:`\Abs{x}`.
* ``\AbsSq``  for the absolute-square, e.g. ``$\AbsSq{\Braket{\Phi}{\Psi}}$`` for :math:`\AbsSq{\Braket{\Phi}{\Psi}}`.
* ``\avg`` for the expectation values, e.g. ``$\avg{\Op{H}}$`` for :math:`\avg{\Op{H}}` (or ``\Avg{}`` for auto-resizing).
* ``\Norm`` for the norm, e.g. ``$\Norm{\ket{\Psi}}$`` for :math:`\Norm{\ket{\Psi}}`.
* ``\identity`` for the identity operator, :math:`\identity`.
* ``\Liouville`` for the Liouvillian symbol, :math:`\Liouville`.
* ``\DynMap`` for the symbolic dynamical map, :math:`\DynMap`.
* ``\dd`` for the differential, e.g. ``$\int f(x) \dd x$`` for :math:`\int f(x) \dd x`.
* Function names / mathematical operators ``\tr``, ``\diag``, ``\abs``, ``\pop``.
* Text labels ``\aux``, ``\opt``, ``\tgt``, ``\init``, ``\lab``, ``\rwa``.

Also see :ref:`math-in-example-notebooks`.

You may use the BibTeX_ plugin for citations.

At any point, from a checkout of the ``krotov`` repository, you may run

.. code-block:: shell

    make docs

or

.. code-block:: shell

   tox -e docs

to generate the documentation locally.

.. _Sphinx: http://www.sphinx-doc.org/en/master/
.. _Restructured Text markup language: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
.. _docstring: https://www.python.org/dev/peps/pep-0257/
.. _"Google Style" format: http://www.sphinx-doc.org/en/master/usage/extensions/example_google.html#example-google
.. _napoleon extension: http://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
.. _mathematical formulas in LaTeX syntax: http://www.sphinx-doc.org/en/1.6/ext/math.html
.. _mathjax: http://www.sphinx-doc.org/en/master/usage/extensions/math.html#module-sphinx.ext.mathjax
.. _BibTeX: https://sphinxcontrib-bibtex.readthedocs.io/en/latest/
.. _Matplotlib Sphinx cheat sheet: https://matplotlib.org/sampledoc/cheatsheet.html


Deploy the documentation
------------------------

The documentation is automatically deployed to
https://qucontrol.github.io/krotov/ (the gh-pages_ associated with the
:mod:`krotov` package's Github repository) every time commits are pushed to
Github. This is done via
`Github Actions <https://github.com/qucontrol/krotov/actions?query=workflow%3ADocs>`_
as configured in the
`workflow file <https://github.com/qucontrol/krotov/blob/master/.github/workflows/docs.yml>`_
at ``.github/workflows/docs.yml``.
The documentation for all versions of :mod:`krotov` is visible on the
``gh-pages`` git branch. Any changes that are committed and pushed from this
branch will be deployed to the online documentation. Do not routinely perform
manual edits on the ``gh-pages`` branch! Let Github Actions do its job of
automatically deploying documentation instead.

.. _gh-pages: https://pages.github.com


.. _ContributeExamples:

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

    .. code-block:: shell

        make notebooks

* Check that the examples can be verified across different Python version by running

    .. code-block:: shell

        make test

* You may also verify that the example is properly integrated in the documentation by running

    .. code-block:: shell

        make docs


.. _math-in-example-notebooks:

Math in Example Notebooks
~~~~~~~~~~~~~~~~~~~~~~~~~

You may use the same tex macros described in the :ref:`write-documentation` section.
However, for the macros to work when viewing the notebook by itself, they must
be redefined locally. To this end, add a markdown cell underneath the top cell
that contains the imported packages (see above). The cell must contain the following:

.. code-block:: tex

    $\newcommand{tr}[0]{\operatorname{tr}}
    \newcommand{diag}[0]{\operatorname{diag}}
    \newcommand{abs}[0]{\operatorname{abs}}
    \newcommand{pop}[0]{\operatorname{pop}}
    \newcommand{aux}[0]{\text{aux}}
    \newcommand{opt}[0]{\text{opt}}
    \newcommand{tgt}[0]{\text{tgt}}
    \newcommand{init}[0]{\text{init}}
    \newcommand{lab}[0]{\text{lab}}
    \newcommand{rwa}[0]{\text{rwa}}
    \newcommand{bra}[1]{\langle#1\vert}
    \newcommand{ket}[1]{\vert#1\rangle}
    \newcommand{Bra}[1]{\left\langle#1\right\vert}
    \newcommand{Ket}[1]{\left\vert#1\right\rangle}
    \newcommand{Braket}[2]{\left\langle #1\vphantom{#2} \mid #2\vphantom{#1}\right\rangle}
    \newcommand{op}[1]{\hat{#1}}
    \newcommand{Op}[1]{\hat{#1}}
    \newcommand{dd}[0]{\,\text{d}}
    \newcommand{Liouville}[0]{\mathcal{L}}
    \newcommand{DynMap}[0]{\mathcal{E}}
    \newcommand{identity}[0]{\mathbf{1}}
    \newcommand{Norm}[1]{\lVert#1\rVert}
    \newcommand{Abs}[1]{\left\vert#1\right\vert}
    \newcommand{avg}[1]{\langle#1\rangle}
    \newcommand{Avg}[1]{\left\langle#1\right\rangle}
    \newcommand{AbsSq}[1]{\left\vert#1\right\vert^2}
    \newcommand{Re}[0]{\operatorname{Re}}
    \newcommand{Im}[0]{\operatorname{Im}}$

Upon executing the cell the definitions will be hidden, but the defined macros
will be available in any cell in the rest of the notebook.

.. _watermark: https://github.com/rasbt/watermark

Versioning
----------

Releases should follow `Semantic Versioning`_, and version numbers published to
PyPI_ must be compatible with :pep:`440`.

In short, versions number follow the pattern `major.minor.patch`, e.g.
``0.1.0`` for the first release, and ``1.0.0`` for the first *stable* release.
If necessary, pre-release versions might be published as e.g:

.. code-block:: none

    1.0.0-dev1  # developer's preview 1 for release 1.0.0
    1.0.0-rc1   # release candidate 1 for 1.0.0

Errors in the release metadata or documentation only may be fixed in a
post-release, e.g.:

.. code-block:: none

    1.0.0.post1  # first post-release after 1.0.0

Post-releases should be used sparingly, but they are acceptable even though
they are not supported by the `Semantic Versioning`_ specification.

The current version is available through the ``__version__`` attribute of the
:mod:`krotov` package:

.. doctest::

    >>> import krotov
    >>> krotov.__version__   # doctest: +SKIP

Between releases, ``__version__`` on the master branch should either be the
version number of the last release, with "+dev" appended (as a
`"local version identifier"`_), or the version number of the next planned
release, with "-dev" appended (`"pre-release identifier"`_ with extra dash).
The version string "1.0.0-dev1+dev" is a valid value after the "1.0.0-dev1"
pre-release. The "+dev" suffix must never be included in a release to PyPI_.

Note that twine_ applies normalization_ to the above recommended forms to
make them strictly compatible with :pep:`440`, before uploading to PyPI_. Users
installing the package through pip_ may use the original version specification
as well as the normalized one (or any other variation that normalizes to the
same result).


Making a Release
----------------

Relesases can only be made by administrators of the Krotov Github repo who are
also listed as Maintainers on https://pypi.org/project/krotov/.

They must have GPG set up to allow for signed commits, and be able to locally
produce documentation artifacts (``make docs-artifacts``).

A release is made by running

.. code-block:: shell

    make release

which executes ``scripts/release.py``. Follow all the prompts.

Releases must be tagged in git, using the version string prefixed by "v",
e.g. ``v1.0.0-dev1`` and ``v1.0.0``. As prompted for by the release script,
after pushing the tag, an official Github-release must be created manually at
https://github.com/qucontrol/krotov/releases, with the proper release notes and
the documentation artifacts as binary attachments.

.. _Semantic Versioning: https://semver.org
.. _"local version identifier": https://www.python.org/dev/peps/pep-0440/#local-version-identifiers
.. _"pre-release identifier": https://www.python.org/dev/peps/pep-0440/#pre-releases
.. _normalization: https://legacy.python.org/dev/peps/pep-0440/#id29
.. _PyPI: http://pypi.org
.. _twine: https://twine.readthedocs.io/en/latest/
.. _pip: https://pip.readthedocs.io/en/stable/


Developers' How-Tos
-------------------

The following assumes your current working directory is a checkout of
``krotov``, and that you have successfully run ``make test`` (which creates
the tox environments that development relies on).


How to run a jupyter notebook server for working on the example notebooks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A notebook server that is isolated to the proper testing environment can be started via the Makefile::

    make jupyter-notebook

This is equivalent to::

    tox -e run-cmd -- jupyter notebook --config=/dev/null

You may run this with your own options, if you prefer. The
``--config=/dev/null`` guarantees that the notebook server is completely
isolated. Otherwise, configuration files from your home directly (see
`Jupyter’s Common Configuration system`_)  may influence the server. Of
course, if you know what you're doing, you may want this.

If you prefer, you may also use the newer jupyterlab::

    make jupyter-lab

How to convert an example notebook to a script for easier debugging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Interactive debugging in notebooks is difficult. It becomes much easier if
you convert the notebook to a script first.  To convert a notebook to an
(I)Python script and run it with automatic debugging, execute e.g.::

    tox -e run-cmd -- jupyter nbconvert --to=python --stdout docs/notebooks/01_example_transmon_xgate.ipynb > debug.py
    tox -e run-cmd -- ipython --pdb debug.py

You can then also set a manual breakpoint by inserting the following line anywhere in the code:

.. code-block:: python

    from IPython.terminal.debugger import set_trace; set_trace() # DEBUG


How to make ``git diff`` work for notebooks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install nbdime_ and run ``nbdime config-git --enable --global`` to `enable the git integration`_.

.. _nbdime: https://nbdime.readthedocs.io/en/latest/index.html
.. _enable the git integration: https://nbdime.readthedocs.io/en/latest/index.html#git-integration-quickstart


How to commit failing tests or example notebooks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The test-suite on the ``master`` branch should always pass without error. If you
would like to commit any example notebooks or tests that currently fail, as a
form of `test-driven development`_, you have two options:

*   Push onto a topic branch (which are allowed to have failing tests), see
    the :ref:`BranchingModel`. The failing tests can then be fixed by
    adding commits to the same branch.

*   Mark the test as failing. For normal tests, add a decorator::

        @pytest.mark.xfail

    See the `pytest documentation on skip and xfail`_ for details.

    For notebooks, the equivalent to the decorator is to add a comment to the
    first line of the failing cell, either::

        # NBVAL_RAISES_EXCEPTION

    (preferably), or::

        # NBVAL_SKIP

    (this may affect subsequent cells, as the marked cell is not executed at all).
    See the `documentation of the nbval pluging on skipping and exceptions`_ for details.


How to run a subset of tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run e.g. only the tests defined in ``tests/test_krotov.py``, use any of the following::

    make test TESTS=tests/test_krotov.py

    tox -e py37-test -- tests/test_krotov.py

    tox -e run-cmd -- pytest tests/test_krotov.py

    .tox/py37/bin/pytest tests/test_krotov.py

See the `pytest test selection docs`_ for details.


How to run only as single test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Decorate the test with e.g. ``@pytest.mark.xxx``, and then run, e.g::

    tox -e run-cmd -- pytest -m xxx tests/

See the `pytest documentation on markers`_ for details.


How to run only the doctests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the following::

    tox -e run-cmd -- pytest --doctest-modules src


How to go into an interactive debugger
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Optionally, install the `pdbpp` package into the testing environment, for a
better experience::

    tox -e run-cmd -- pip install pdbpp

Then:

- before the line where you went to enter the debugger, insert a line::

    from IPython.terminal.debugger import set_trace; set_trace() # DEBUG

- Run ``pytest`` with the option ``-s``, e.g.::

    tox -e run-cmd -- pytest -m xxx -s tests/

You may also see the `pytest documentation on automatic debugging`_.


How to see the debug logger output in the example notebooks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :func:`.optimize_pulses` routine generates some logger messages for
debugging purposes. To see these messages, set the level of "krotov" logger to
INFO or DEBUG:

.. code-block:: python

   import logging
   logger = logging.getLogger('krotov')
   logger.setLevel(logging.DEBUG)


You can also configure the logger with custom formatters, e.g. to show the
messages with time stamps:

.. code-block:: python

   ch = logging.StreamHandler()
   ch.setLevel(logging.INFO)
   formatter = logging.Formatter("%(asctime)s:%(message)s")
   ch.setFormatter(formatter)
   logger.addHandler(ch)
   logging.getLogger().handlers = [] # disable root handlers


See the `Configure Logging`_ section of the Python documentation for more details.


How to use quantum mechanical tex macros
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For docstrings or ``*.rst`` files, see :ref:`write-documentation`. For notebooks, see :ref:`math-in-example-notebooks`.


.. _Jupyter’s Common Configuration system: https://jupyter-notebook.readthedocs.io/en/stable/config_overview.html#jupyter-s-common-configuration-system
.. _Closing issues using keywords: https://help.github.com/articles/closing-issues-using-keywords/
.. _pytest test selection docs: https://docs.pytest.org/en/latest/usage.html#specifying-tests-selecting-tests
.. _pytest documentation on markers: https://docs.pytest.org/en/latest/example/markers.html
.. _pytest documentation on automatic debugging: https://docs.pytest.org/en/latest/usage.html#dropping-to-pdb-python-debugger-on-failures
.. _test-driven development: https://en.wikipedia.org/wiki/Test-driven_development
.. _pytest documentation on skip and xfail: https://docs.pytest.org/en/latest/skipping.html
.. _documentation of the nbval pluging on skipping and exceptions: https://nbval.readthedocs.io/en/latest/#Skipping-specific-cells
.. _Configure Logging: https://docs.python.org/3/howto/logging.html#configuring-logging
