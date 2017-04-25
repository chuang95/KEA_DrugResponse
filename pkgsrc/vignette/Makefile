build:
	R CMD Rscript -e "library(knitr); knit('KEADrugResponse.Rnw')"
	latex KEADrugResponse.tex
	bibtex KEADrugResponse.aux
	latex KEADrugResponse.tex
	latex KEADrugResponse.tex
	dvipdf KEADrugResponse.dvi

clean:
	rm -rf *.tex *.dvi *.pdf

