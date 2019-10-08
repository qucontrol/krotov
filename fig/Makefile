PDFLATEXOPTS = -file-line-error -interaction=nonstopmode -halt-on-error
ALL = krotovscheme.pdf oct_decision_tree.pdf tls_state_to_state_iters.pdf
PYTHON = python

all: $(ALL)

%.pdf: %.tex
	@echo "\n** generate figure: $< -> $@"
	@pdflatex ${PDFLATEXOPTS} $<

%.pdf: %.py
	@echo "\n** generate figure: $< -> $@"
	@$(PYTHON) $<

%.tex: %.py
	@$(PYTHON) $<

tls_state_to_state_iters.pdf: tls_state_to_state_iters.py mpl.py matplotlibrc

krotovscheme.pdf: krotovscheme.tex

oct_decision_tree.pdf: oct_decision_tree.tex

clean:
	@rm -f *.aux
	@rm -f *.log
	@rm -f *.toc
	@rm -f *.blg
	@rm -rf *.out
	@rm -f *.bak
	@rm -f *.ilg
	@rm -f *.snm
	@rm -f *.nav
	@rm -f *.table
	@rm -f *.dvi
	@rm -f *.fls
	@rm -f *.bbl
	@rm -f *~
	@rm -f *Notes.bib
	@rm -f *-eps-converted-to.pdf
	@rm -f *.fdb_latexmk
	@rm -f *.synctex.gz*
	@rm -f ./latexrun

distclean:
	rm -f $(ALL)

.PHONY: all clean distclean
