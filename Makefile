all: distance-estimator.html

clean:
	rm -f distance-estimator.html

.PHONY: all clean

# Render Markdown to HTML
%.html: %.md
	Rscript -e 'rmarkdown::render("$<", output_format = "html_document")'

# Render Markdown to PDF
%.pdf: %.md
	Rscript -e 'rmarkdown::render("$<", output_format = "pdf_document")'

# Render LaTeX to PDF
%.pdf: %.tex
	pdflatex -halt-on-error $<
