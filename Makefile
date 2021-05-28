.PHONY: black black-check clean clean-build clean-tests clean-venv coverage dist dist-check docs docs-pdf help install isort isort-check jupyter-lab jupyter-notebook flake8-check pylint-check notebooks pre-commit-hooks release spellcheck test test-upload uninstall upload
.DEFAULT_GOAL := help
TOXOPTIONS ?=
TOXINI ?= tox.ini
TOX = tox -c $(TOXINI) $(TOXOPTIONS)
# and empty TESTS delegates to tox.ini
TESTS ?=


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


help:  ## show this help
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)


bootstrap: ## verify that tox is available and pre-commit hooks are active
	python scripts/bootstrap.py

clean: ## remove all build, test, coverage, and Python artifacts, as well as environments
	python scripts/clean.py all

clean-build: ## remove build artifacts
	python scripts/clean.py build

clean-tests: ## remove test and coverage artifacts
	python scripts/clean.py tests

clean-venv: ## remove tox virtual environments
	python scripts/clean.py venv

clean-docs: ## remove documentation artifacts
	python scripts/clean.py docs

flake8-check: bootstrap ## check style with flake8
	$(TOX) -e run-flake8

pylint-check: bootstrap ## check style with pylint
	$(TOX) -e run-pylint

test: bootstrap ## run tests for current stable Python release
	$(TOX) -e py38-test -- $(TESTS)

test36: bootstrap ## run tests for Python 3.6
	$(TOX) -e py36-test -- $(TESTS)

test37: bootstrap ## run tests for Python 3.7
	$(TOX) -e py37-test -- $(TESTS)

test38: bootstrap ## run tests for Python 3.8
	$(TOX) -e py38-test -- $(TESTS)

docs: bootstrap ## generate Sphinx HTML documentation, including API docs
	$(TOX) -e docs
	@echo "open docs/_build/html/index.html"

docs-pdf: bootstrap  ## generate a PDF version of the documentation
	$(TOX) -e docs -- _build/tex -b latex -d _build/doctree
	cp docs/krotovscheme.pdf docs/oct_decision_tree.pdf docs/_build/tex/
	$(TOX) -e run-cmd -- python docs/build_pdf.py

KROTOV_VERSION = $(shell grep __version__ src/krotov/__init__.py | sed "s/__version__ = '//g" | sed "s/'//g")

docs-artifacts: docs docs-pdf  ## create the documentation artifacts for a release
	mkdir -p docs/_build/artifacts
	$(TOX) -e run-cmd -- zip-folder docs/_build/html --root-folder=krotov-v$(KROTOV_VERSION) --outfile docs/_build/artifacts/krotov-v$(KROTOV_VERSION).zip
	cp docs/_build/tex/krotov.pdf docs/_build/artifacts/krotov-v$(KROTOV_VERSION).pdf


spellcheck: bootstrap ## check spelling in docs
	$(TOX) -e run-cmd -- pip install sphinxcontrib-spelling
	SPELLCHECK=en_US $(TOX) -e docs -- _build/spell -b spelling -d _build/doctree

black-check: bootstrap ## Check all src and test files for complience to "black" code style
	$(TOX) -e run-blackcheck

black: bootstrap ## Apply 'black' code style to all src and test files
	$(TOX) -e run-black

isort-check: bootstrap ## Check all src and test files for correctly sorted imports
	$(TOX) -e run-isortcheck

isort: bootstrap ## Sort imports in all src and test files
	$(TOX) -e run-isort

coverage: test38  ## generate coverage report in ./htmlcov
	$(TOX) -e coverage
	@echo "open htmlcov/index.html"

test-upload: bootstrap clean-build dist ## package and upload a release to test.pypi.org
	$(TOX) -e run-cmd -- twine check dist/*
	$(TOX) -e run-cmd -- twine upload --repository-url https://test.pypi.org/legacy/ dist/*

upload: bootstrap clean-build dist ## package and upload a release to pypi.org
	$(TOX) -e run-cmd -- twine check dist/*
	$(TOX) -e run-cmd -- twine upload dist/*

release: bootstrap ## Create a new version, package and upload it
	python3.7 -m venv .venv/release
	.venv/release/bin/python -m pip install click gitpython pytest
	.venv/release/bin/python ./scripts/release.py

dist: bootstrap ## builds source and wheel package
	$(TOX) -e run-cmd -- python setup.py sdist
	$(TOX) -e run-cmd -- python setup.py bdist_wheel
	ls -l dist

dist-check: bootstrap ## Check all dist files for correctness
	$(TOX) -e run-cmd -- twine check dist/*

install: clean-build ## install the package to the active Python's site-packages
	pip install .

uninstall:  ## uninstall the package from the active Python's site-packages
	pip uninstall krotov

# How to execute notebook files
%.ipynb.log: %.ipynb
	@echo ""
	$(TOX) -e run-cmd -- jupyter nbconvert --to notebook --execute --inplace --allow-errors --ExecutePreprocessor.kernel_name='python3'  --ExecutePreprocessor.timeout=-1 --config=/dev/null $< 2>&1 | tee $@

NOTEBOOKFILES = $(shell find docs/notebooks -maxdepth 1 -iname '*.ipynb')
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
