.PHONY: clean clean-test clean-pyc clean-build clean-venvs line pep8 docs dist install develop help
.DEFAULT_GOAL := help
CONDA_PACKAGES = qutip
TESTENV =
#TESTENV = MATPLOTLIBRC=tests
TESTOPTIONS = --doctest-modules --cov=krotov --nbval --sanitize-with docs/nbval_sanitize.cfg
TESTS = src tests docs/notebooks/*.ipynb docs/*.rst docs/*.ipynb


define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test clean-venvs ## remove all build, test, coverage, and Python artifacts, as well as environments
	$(MAKE) -C docs clean

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	rm -fr src/krotov.egg-info
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

lint: ## check style with flake8
	flake8 src tests

pep8: ## check style with pep8
	pep8 src tests


test:  test35 test36 ## run tests on every Python version


.venv/py35/bin/py.test:
	@conda create -y -m --override-channels -c defaults -p .venv/py35 python=3.5
	@# if the conda installation does not work, simply comment out the following line, and let pip handle it
	@conda install -y --override-channels -c defaults -c conda-forge -p .venv/py35 $(CONDA_PACKAGES)
	@.venv/py35/bin/pip install -e .[dev]

test35: .venv/py35/bin/py.test ## run tests for Python 3.5
	$(TESTENV) $< -v $(TESTOPTIONS) $(TESTS)



.venv/py36/bin/py.test:
	@conda create -y -m --override-channels -c defaults -p .venv/py36 python=3.6
	@# if the conda installation does not work, simply comment out the following line, and let pip handle it
	@conda install -y --override-channels -c defaults -c conda-forge -p .venv/py36 $(CONDA_PACKAGES)
	@.venv/py36/bin/pip install -e .[dev]


test36: .venv/py36/bin/py.test ## run tests for Python 3.6
	$(TESTENV) $< -v $(TESTOPTIONS) $(TESTS)


.venv/py36/bin/python: .venv/py36/bin/py.test

.venv/py36/bin/sphinx-build: .venv/py36/bin/py.test

.venv/py36/bin/jupyter: .venv/py36/bin/py.test

# How to execute notebook files
%.ipynb.log: %.ipynb .venv/py36/bin/jupyter
	@echo ""
	@.venv/py36/bin/jupyter nbconvert --to notebook --execute --inplace --allow-errors --ExecutePreprocessor.timeout=180 --ExecutePreprocessor.kernel_name='python3' --config=/dev/null $< 2>&1 | tee $@

NOTEBOOKFILES = $(shell find docs/notebooks/ -iname '*.ipynb'  -maxdepth 1) docs/10_time_discretization.ipynb
NOTEBOOKLOGS = $(patsubst %.ipynb,%.ipynb.log,$(NOTEBOOKFILES))

notebooks: $(NOTEBOOKLOGS)  ## re-evaluate the notebooks in docs/notebooks
	@echo ""
	@echo "All notebook are now up to date; the were executed using the python3 kernel"
	@.venv/py36/bin/jupyter kernelspec list | grep python3


docs: .venv/py36/bin/sphinx-build ## generate Sphinx HTML documentation, including API docs
	$(MAKE) -C docs SPHINXBUILD=../.venv/py36/bin/sphinx-build clean
	$(MAKE) -C docs SPHINXBUILD=../.venv/py36/bin/sphinx-build html
	@echo "open docs/_build/html/index.html"

spellcheck: .venv/py36/bin/sphinx-build ## check spelling in docs
	@.venv/py36/bin/pip install sphinxcontrib-spelling
	SPELLCHECK=en_US $(MAKE) -C docs SPHINXBUILD=../.venv/py36/bin/sphinx-build spelling

coverage: test36  ## generate coverage report in ./htmlcov
	.venv/py36/bin/coverage html
	@echo "open htmlcov/index.html"

test-upload: .venv/py36/bin/python clean-build clean-pyc dist ## package and upload a release to test.pypi.org
	.venv/py36/bin/twine check dist/*
	.venv/py36/bin/twine upload dist/*

upload: .venv/py36/bin/python clean-build clean-pyc dist ## package and upload a release to pypi.org
	.venv/py36/bin/twine check dist/*
	.venv/py36/bin/twine upload --repository-url https://test.pypi.org/legacy/ dist/*

release: clean .venv/py36/bin/python ## Create a new version, package and upload it
	.venv/py36/bin/python ./scripts/release.py krotov


dist: .venv/py36/bin/python clean-build clean-pyc ## builds source and wheel package
	@$< setup.py sdist
	@$< setup.py bdist_wheel
	ls -l dist

dist-check: .venv/py36/bin/python  ## Check all dist files for correctness
	.venv/py36/bin/twine check dist/*

install: clean-build clean-pyc ## install the package to the active Python's site-packages
	pip install .

uninstall:  ## uninstall the package from the active Python's site-packages
	pip uninstall krotov

develop: clean-build clean-pyc ## install the package to the active Python's site-packages, in develop mode
	pip install -e .

develop-test: develop ## run tests within the active Python environment
	$(TESTENV) py.test -v $(TESTOPTIONS) $(TESTS)

develop-docs: develop  ## generate Sphinx HTML documentation, including API docs, within the active Python environment
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	@echo "open docs/_build/html/index.html"

jupyter-notebook: .venv/py36/bin/jupyter  ## run a notebook server for editing the examples
	.venv/py36/bin/jupyter notebook --config=/dev/null

jupyter-lab: .venv/py36/bin/jupyter  ## run a jupyterlab server for editing the examples
	@.venv/py36/bin/pip install -e .[extras]
	.venv/py36/bin/jupyter lab --config=/dev/null
