.PHONY: help init docs docs-artifacts serve clean clean-build clean-tests clean-venv clean-docs distclean test coverage black-check black flake8-check pylint-check listenvs shell upload test-upload release release-test dist dist-check install uninstall notebooks jupyter-lab

TESTS ?= src tests docs

help: init  ## Show this help
	@grep -E '^([a-zA-Z_-]+):.*## ' $(MAKEFILE_LIST) | awk -F ':.*## ' '{printf "%-20s %s\n", $$1, $$2}'

init:  .git/hooks/pre-commit  ## Initialize the repository after cloning

.git/hooks/pre-commit:
	@pre-commit install

docs:  ## Build the documentation
	hatch run docs:build

serve: docs  ## Build and serve the documentation
	hatch run docs:serve

clean: ## Clean up build/doc/testing artifacts
	python scripts/clean.py all

clean-build: ## remove build artifacts
	python scripts/clean.py build

clean-tests: ## remove test and coverage artifacts
	python scripts/clean.py tests

clean-venv: ## remove tox virtual environments
	python scripts/clean.py venv

clean-docs: ## remove documentation artifacts
	python scripts/clean.py docs

distclean: clean  ## Restore to a clean checkout state
	rm -f .git/hooks/pre-commit
	hatch env prune

test: init ## Run the tests (with coverage)
	hatch run cov -- $(TESTS)

coverage: ## Generate coverage report in ./htmlcov (after `make test`)
	hatch run cov-html
	@echo "open htmlcov/index.html"

black-check: ## Check all src and test files for complience to "black" code style
	hatch run lint:black-check

black: ## Apply 'black' code style to all src and test files
	hatch run lint:black

isort-check: ## Check all src and test files for correctly sorted imports
	hatch run lint:isort-check

isort: ## Sort imports in all src and test files
	hatch run lint:isort

flake8-check: ## check style with flake8
	hatch run lint:flake8

pylint-check: ## check style with pylint
	hatch run lint:pylint

listenvs: ## Show all hatch environments
	hatch env show

shell: ## Open a shell with all development tools
	hatch shell

upload: dist ## Package and upload a release to pypi.org
	hatch publish

test-upload: dist ## Package and upload a release to test.pypi.org
	hatch publish --repo test

release: ## Create a new version, package and upload it
	hatch -e release run -- python ./scripts/release.py

release-test: ## Run the release script is "test" mode
	hatch -e release run -- python ./scripts/release.py --test

dist: ## Build source and wheel package
	hatch build --clean

dist-check: ## Check all dist files for correctness
	hatch run release:twine check dist/*

install: clean-build ## install the package to the active Python's site-packages
	pip install .

uninstall:  ## uninstall the package from the active Python's site-packages
	pip uninstall krotov


# How to execute notebook files
%.ipynb.log: %.ipynb
	@echo ""
	hatch run -- jupyter nbconvert --to notebook --execute --inplace --allow-errors --ExecutePreprocessor.kernel_name='python3'  --ExecutePreprocessor.timeout=-1 --config=/dev/null $< 2>&1 | tee $@

NOTEBOOKFILES = $(shell find docs/notebooks -maxdepth 1 -iname '*.ipynb')
NOTEBOOKLOGS = $(patsubst %.ipynb,%.ipynb.log,$(NOTEBOOKFILES))

notebooks: $(NOTEBOOKLOGS)  ## re-evaluate the notebooks
	@echo ""
	@echo "All notebook are now up to date; the were executed using the python3 kernel"
	hatch run -- jupyter kernelspec list | grep python3

jupyter-lab:  ## Run a Jupyterlab server for editing the examples
	hatch run -- jupyter lab
