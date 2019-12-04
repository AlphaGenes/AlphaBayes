#!/bin/sh

TMPS="aux ilg ind log maf mtc* stc* tex toc"

for TMP in $TMPS; do
  rm -f AlphaBayes.${TMP}
done

rst2latex.py --stylesheet preamble.sty AlphaBayes.txt AlphaBayes.tex

pdflatex AlphaBayes.tex
bibtex AlphaBayes.tex
makeindex AlphaBayes.tex
pdflatex AlphaBayes.tex
pdflatex AlphaBayes.tex

for TMP in $TMPS; do
  rm -f AlphaBayes.${TMP}
done

open AlphaBayes.pdf
