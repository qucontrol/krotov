# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = src/krotov
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile


# the krotovscheme images only get re-generated when they are deleted. We don't
# really want ot add a dev dependency on latex and Poppler

krotovscheme.svg:
	pdflatex krotovscheme.tex
	pdftocairo -svg -expand -paperh 675 -paperw 272  krotovscheme.pdf

krotovscheme.pdf:
	pdflatex krotovscheme.tex


clean:
	@-rm -rf $(BUILDDIR)/*
	@-rm -rf API/*.rst
	@-rm -f notebooks/*.log
	@-rm -f *.log
	@-rm -f readme.rst
	@-rm -f krotovscheme.aux krotovscheme.log

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@find . -type d -iname .ipynb_checkpoints -exec rm -rf "{}" \; | true
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
