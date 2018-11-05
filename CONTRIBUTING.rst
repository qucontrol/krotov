.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/qucontrol/krotov/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs / Implement Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for bugs or feature requests. Anybody is welcome to submit a pull request for open issues.


Write Documentation
~~~~~~~~~~~~~~~~~~~

Krotov could always use more documentation, whether as part of the
official Krotov docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/qucontrol/krotov/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Follow `Aaron Meurer's Git Workflow Notes`_ (with ``qucontrol/krotov`` instead of ``sympy/sympy``)

In short, *if you are not a member of the `qucontrol organization`_*,

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

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. Check https://travis-ci.org/qucontrol/krotov/pull_requests
   and make sure that the tests pass for all supported Python versions.

