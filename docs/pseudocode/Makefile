PDFLATEXOPTS = -file-line-error -interaction=nonstopmode -halt-on-error -synctex=1
NAME = krotov_pseudocode

all: ${NAME}.pdf

.PHONY: FORCE
${NAME}.pdf: FORCE ${NAME}.tex
	@echo ""
	@echo "*********"
	@echo "Compiling Main File with ..."
	pdflatex ${PDFLATEXOPTS} ${NAME}.tex
	pdflatex ${PDFLATEXOPTS} ${NAME}.tex
	@echo "Done"


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
	@echo "Done"


.PHONY: all clean
