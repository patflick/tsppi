all: reqs.pdf


reqs.pdf: reqs.tex biblio.bib
	pdflatex reqs.tex
	bibtex reqs.aux
	pdflatex reqs.tex
	pdflatex reqs.tex
