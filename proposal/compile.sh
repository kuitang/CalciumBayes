#/bin/sh

PATH=/usr/local/texlive/2010/bin/x86_64-darwin:$PATH
pdflatex main
bibtex main
pdflatex main
pdflatex main
