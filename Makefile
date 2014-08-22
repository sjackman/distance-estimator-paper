all: mle.pdf

clean:
	rm -f mle.pdf

.PHONY: all clean

%.pdf: %.tex
	pdflatex -halt-on-error $<
