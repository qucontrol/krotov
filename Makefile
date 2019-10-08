PDFLATEXOPTS = -file-line-error -interaction=nonstopmode -halt-on-error -synctex=1
NAME = krotov
FIGS = fig1.pdf fig2.pdf fig3.pdf
EXAMPLES = examples/tls_state_to_state.py examples/quantum_gate_hilbert.py examples/quantum_gate_liouville.py  examples/ensemble.py examples/second_order.py
EXAMPLES_OUT = examples/tls_state_to_state.out examples/quantum_gate_hilbert.out examples/quantum_gate_liouville.out  examples/ensemble.out examples/second_order.out
ARXIVEXTRA = SciPost.cls
PYTHON = .venv/bin/python

all: ${NAME}.pdf

.PHONY: FORCE
${NAME}.pdf: FORCE ${NAME}.tex $(FIGS) $(EXAMPLES) $(EXAMPLES_OUT)
	@echo ""
	@echo "*********"
	@echo "Compiling Main File with ..."
	pdflatex ${PDFLATEXOPTS} ${NAME}.tex
	bibtex ${NAME}
	pdflatex ${PDFLATEXOPTS} ${NAME}.tex
	pdflatex ${PDFLATEXOPTS} ${NAME}.tex
	@echo "Done"


${NAME}.bbl: ${NAME}.pdf

fig/%.pdf: $(PYTHON)
	$(MAKE) -C fig PYTHON=../$(PYTHON)

fig/krotovscheme.pdf: fig/krotovscheme.tex

fig/oct_decision_tree.pdf: fig/oct_decision_tree.tex

fig/tls_state_to_state_iters.pdf: fig/tls_state_to_state_iters.py fig/mpl.py fig/matplotlibrc examples/tls_state_to_state.dump

fig1.pdf: fig/tls_state_to_state_iters.pdf
	cp $< $@

fig2.pdf: fig/krotovscheme.pdf
	cp $< $@

fig3.pdf: fig/oct_decision_tree.pdf
	cp $< $@

$(PYTHON):
	python3.7 -m venv .venv
	$(PYTHON) -m pip install krotov[dev]==1.0.0

%.out: %.py $(PYTHON)
	$(PYTHON) $< | tee $@

test: $(PYTHON)
	$(PYTHON) -m pytest -s --verbose examples/test_examples.py

examples/tls_state_to_state.out: examples/tls_state_to_state.py

examples/quantum_gate_hilbert.out: examples/quantum_gate_hilbert.py

examples/quantum_gate_liouville.out: examples/quantum_gate_liouville.py

examples/ensemble.out: examples/ensemble.py

examples/second_order.out: examples/second_order.py

arxiv.tgz:  ${NAME}.tex ${NAME}.bbl ${FIGS} ${EXAMPLES} ${EXAMPLES_OUT} ${ARXIVEXTRA}
	@rm -rf arxiv
	@mkdir -p arxiv
	@cp ${NAME}.tex ${NAME}.bbl ${FIGS} ${ARXIVEXTRA} arxiv/
	@cp -r examples arxiv/
	@COPYFILE_DISABLE=1 tar -czf arxiv.tgz -C arxiv .

clean:
	@echo "Cleaning up files from LaTeX compilation ..."
	@rm -f *.aux
	@rm -f *.log
	@rm -f *.toc
	@rm -f *.blg
	@rm -rf *.out
	@rm -rf .cache
	@rm -rf examples/__pycache__
	@rm -f *.bak
	@rm -f *.ilg
	@rm -f *.snm
	@rm -f *.nav
	@rm -f *.table
	@rm -f *.dvi
	@rm -f *.fls
	@rm -f *~
	@rm -f *Notes.bib
	@rm -f *-eps-converted-to.pdf
	@rm -f *.fdb_latexmk
	@rm -f *.synctex.gz*
	@rm -f .latexrun.db*
	make -C fig clean
	@echo "Done"

distclean: clean
	@echo "Removing all compiled files"
	@rm -f ${NAME}.pdf
	@rm -f *.bbl
	@rm -f $(FIGS)
	@rm -f latexrun
	@make -C fig distclean
	@rm -rf submit
	@rm -rf arxiv
	@rm -f submit.tgz
	@rm -f arxiv.tgz
	@rm -rf .venv
	@echo "Done"

.PHONY: all update clean distclean
