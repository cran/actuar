### -*-Makefile-*- to build actuar vignettes
##
##  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

## List of vignette to build
VIGNETTES = actuar.pdf coverage.pdf credibility.pdf \
	distributions.pdf modeling.pdf risk.pdf \
	simulation.pdf

## Toolset
SWEAVE = "$(R_HOME)/bin/R" CMD Sweave --encoding="utf-8"
TEXI2DVI = LATEX=xelatex TEXINDY=makeindex texi2dvi -b
RM = rm -rf

all: pdf

%.pdf: %.tex
	${TEXI2DVI} '$<'

.PHONY: pdf
pdf: ${VIGNETTES}

.PHONY: clean
clean:
	${RM} *.tex *-[0-9][0-9][0-9].pdf \
	      *.aux *.bbl *.blg *.log *.out *~ Rplots* \
	      auto/
