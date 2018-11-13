# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = src/krotov
SOURCEDIR     = .
BUILDDIR      = _build

krotovscheme.svg: krotovscheme.pdf
	pdftocairo -svg -expand -paperh 675 -paperw 272  krotovscheme.pdf

krotovscheme.pdf:
	pdflatex krotovscheme.tex


# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

clean:
	@-rm -rf $(BUILDDIR)/*
	@-rm -rf API/*.rst
	@-rm -f notebooks/*.log
	@rm -f krotovscheme.aux krotovscheme.log

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
