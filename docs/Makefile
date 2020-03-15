# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = src/krotov
SOURCEDIR     = .
BUILDDIR      = _build/html

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile figs


# The krotovscheme images only get re-generated when they are deleted. We don't
# really want ot add a dev dependency on latex and Poppler
#
# Note: the figures are in landscape mode. Thus, paperw and paperh are
# switched. We also apply a scaling figure, to get a bigger svg figure.

figs: krotovscheme.pdf krotovscheme.svg oct_decision_tree.pdf oct_decision_tree.svg

krotovscheme.svg:
	pdflatex krotovscheme.tex
	pdftocairo -svg -expand -paperw 340 -paperh 674  krotovscheme.pdf  # scaling factor 2

krotovscheme.pdf:
	pdflatex krotovscheme.tex


oct_decision_tree.svg:
	pdflatex oct_decision_tree.tex
	pdftocairo -svg -expand -paperw 400 -paperh 744  oct_decision_tree.pdf  # scaling factor 1.5

oct_decision_tree.pdf:
	pdflatex oct_decision_tree.tex


clean:
	@python clean.py

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@find . -type d -iname .ipynb_checkpoints -exec rm -rf "{}" \; | true
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
