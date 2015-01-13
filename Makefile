all: distance-estimator.html distance-estimator.pdf

clean:
	rm -f \
		distance-estimator.html \
		distance-estimator.pdf \
		distance-estimator.tex

install-deps:
	brew install pandoc
	brew cask install mactex

.PHONY: all clean
.DELETE_ON_ERROR:
.SECONDARY:

mathjax=https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML

# Render Markdown to HTML
%.html: %-before.md %.md
	pandoc --template=template -F pandoc-citeproc --mathjax=$(mathjax) -o $@ $^

# Render Markdown to LaTeX
%.tex: %.md
	pandoc --template=template -F pandoc-citeproc -o $@ $<

# Render LaTeX to PDF
%.pdf: %.tex
	pdflatex -halt-on-error $<

# Download peerj.csl
peerj.csl:
	curl -O https://raw.githubusercontent.com/citation-style-language/styles/master/peerj.csl
