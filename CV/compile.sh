rm -f CV.{aux,bbl,blg,log,out} &&
pdflatex CV.tex &&
bibtex CV && 
pdflatex CV.tex &&
pdflatex CV.tex 

