all: pdf 

pdf: DClusterm.pdf

%.pdf: %.tex
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*
	pdflatex $*

clean:
	rm *tex *aux *log *bbl *blg *out jss_paper-*.pdf jss.* 
	
