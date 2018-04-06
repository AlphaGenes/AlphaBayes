#!/bin/sh

rm -f \
AlphaBayes.aux \
AlphaBayes.ilg \
AlphaBayes.ind \
AlphaBayes.log \
AlphaBayes.maf \
AlphaBayes.mtc \
AlphaBayes.mtc0 \
AlphaBayes.pdf \
AlphaBayes.stc1 \
AlphaBayes.stc2 \
AlphaBayes.stc3 \
AlphaBayes.stc4 \
AlphaBayes.stc5 \
AlphaBayes.stc6 \
AlphaBayes.stc7 \
AlphaBayes.stc8 \
AlphaBayes.stc9 \
AlphaBayes.tex \
AlphaBayes.toc

rst2latex.py --stylesheet preamble.sty AlphaBayes.txt AlphaBayes.tex

pdflatex AlphaBayes.tex
bibtex AlphaBayes.tex
makeindex AlphaBayes.tex
pdflatex AlphaBayes.tex
pdflatex AlphaBayes.tex

open AlphaBayes.pdf
