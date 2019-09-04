.PHONY: black black-check clean clean-build clean-pyc clean-test clean-venvs coverage dist dist-check docs help install isort isort-check jupyter-lab jupyter-notebook flake8-check pylint-check notebooks pre-commit-hooks release spellcheck test test-upload uninstall upload
.DEFAULT_GOAL := help
TOXOPTIONS =
TOXINI = tox.ini
TOX = tox -c $(TOXINI) $(TOXOPTIONS)


define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
    match = re.match(r'^([a-z0-9A-Z_-]+):.*?## (.*)$$', line)
    if match:
        target, help = match.groups()
        print("%-20s %s" % (target, help))
print("""
All make commands delegate to tox. Environments will be created with venv by
default. Alternatively, conda can be used by using the tox-conda.ini file,
e.g. with `make TOXINI=tox-conda.ini test` You may also run `tox` directly. See
`tox -av` for a list of available environments.
""")
endef
export PRINT_HELP_PYSCRIPT

define BOOTSTRAP_PYSCRIPT
import sys
import os
from textwrap import dedent
try:
    import tox
    if not os.path.isfile(".git/hooks/pre-commit"):
        print("bootstrapping pre-commit hook")
        tox.cmdline(['-e', 'run-cmd', '--', 'pre-commit', 'install'])
except ImportError:
    print(dedent("""
    tox is not available. See https://tox.readthedocs.io for installation
    instructions.
    """))
    sys.exit(1)
endef
export BOOTSTRAP_PYSCRIPT

help:  ## show this help
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)


bootstrap: ## verify that tox is available and pre-commit hooks are active
	@python -c "$$BOOTSTRAP_PYSCRIPT"

clean: clean-docs clean-build clean-pyc clean-test clean-venvs ## remove all build, test, coverage, and Python artifacts, as well as environments

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	rm -fr src/*.egg-info
	rm -fr pip-wheel-metadata
	find tests src -name '*.egg-info' -exec rm -fr {} +
	find tests src -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find tests src -name '*.pyc' -exec rm -f {} +
	find tests src -name '*.pyo' -exec rm -f {} +
	find tests src -name '*~' -exec rm -f {} +
	find tests src -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -f .coverage
	rm -fr htmlcov/

clean-venvs: ## remove testing/build environments
	rm -fr .tox
	rm -fr .venv

clean-docs: ## remove documentation artifacts
	$(MAKE) -C docs clean

flake8-check: ## check style with flake8
	$(TOX) -e run-flake8

pylint-check: ## check style with pylint
	$(TOX) -e run-pylint

test: bootstrap ## run tests on every supported Python version
	$(TOX) -e py35-test,py36-test,py37-test

test35: bootstrap ## run tests for Python 3.5
	$(TOX) -e py35-test

test36: bootstrap ## run tests for Python 3.6
	$(TOX) -e py36-test

test37: bootstrap ## run tests for Python 3.7
	$(TOX) -e py37-test

docs: bootstrap ## generate Sphinx HTML documentation, including API docs
	$(TOX) -e docs
	@echo "open docs/_build/index.html"

spellcheck: bootstrap ## check spelling in docs
	$(TOX) -e run-cmd -- pip install sphinxcontrib-spelling
	SPELLCHECK=en_US $(TOX) -e docs -- -b spelling

black-check: bootstrap ## Check all src and test files for complience to "black" code style
	$(TOX) -e run-blackcheck

black: bootstrap ## Apply 'black' code style to all src and test files
	$(TOX) -e run-black

isort-check: bootstrap ## Check all src and test files for correctly sorted imports
	$(TOX) -e run-isortcheck

isort: bootstrap ## Sort imports in all src and test files
	$(TOX) -e run-isort

coverage: test37  ## generate coverage report in ./htmlcov
	$(TOX) -e coverage
	@echo "open htmlcov/index.html"

test-upload: bootstrap clean-build clean-pyc dist ## package and upload a release to test.pypi.org
	$(TOX) -e run-cmd -- twine check dist/*
	$(TOX) -e run-cmd -- twine upload --repository-url https://test.pypi.org/legacy/ dist/*

upload: bootstrap clean-build clean-pyc dist ## package and upload a release to pypi.org
	$(TOX) -e run-cmd -- twine check dist/*
	$(TOX) -e run-cmd -- twine upload dist/*

release: bootstrap ## Create a new version, package and upload it
	$(TOX) -e run-cmd -- python ./scripts/release.py

dist: bootstrap ## builds source and wheel package
	$(TOX) -e run-cmd -- python setup.py sdist
	$(TOX) -e run-cmd -- python setup.py bdist_wheel
	ls -l dist

dist-check: bootstrap ## Check all dist files for correctness
	$(TOX) -e run-cmd -- twine check dist/*

install: clean-build clean-pyc ## install the package to the active Python's site-packages
	pip install .

uninstall:  ## uninstall the package from the active Python's site-packages
	pip uninstall krotov

# How to execute notebook files
%.ipynb.log: %.ipynb
	@echo ""
	$(TOX) -e run-cmd -- jupyter nbconvert --to notebook --execute --inplace --allow-errors --ExecutePreprocessor.kernel_name='python3' --config=/dev/null $< 2>&1 | tee $@

NOTEBOOKFILES = $(shell find docs/ -maxdepth 1 -iname '*.ipynb')
NOTEBOOKLOGS = $(patsubst %.ipynb,%.ipynb.log,$(NOTEBOOKFILES))

notebooks: bootstrap $(NOTEBOOKLOGS)  ## re-evaluate the notebooks
	@echo ""
	@echo "All notebook are now up to date; the were executed using the python3 kernel"
	$(TOX) -e run-cmd -- jupyter kernelspec list | grep python3

jupyter-notebook: bootstrap ## run a notebook server for editing the examples
	$(TOX) -e run-cmd -- jupyter notebook --config=/dev/null

jupyter-lab: bootstrap ## run a jupyterlab server for editing the examples
	$(TOX) -e run-cmd -- pip install jupyterlab
	$(TOX) -e run-cmd -- jupyter lab --config=/dev/null
