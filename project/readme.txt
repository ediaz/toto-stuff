This folder contains the Python code for the 
2D variable velcity Helmhotz solver.

I tried to make my report fully reproducible, 
so if interest you can create the figures
compile the tex files.

Requirements:
matplotlib
scipy
numpy
pdflatex

Content:
This folder has 5 main Python codes:
1- plot.py:
   contains plotting wrappers based on matplotlib

2- difference.py:
  contain a class with the implementation of Fornberg, 1988 paper
  for digital differentiators

3- helmhotzLibrary.py:
  contains the class for the Helmhotz solver (class helmhotz)
  and the class for the PML weights and its derivatives.

4- rsflikeLibrary:
  defines objects that have an array of floats and grid 
  geometry atributes. I tried to mimic the Madagascar 
  format 'rsf'.

5- helmhotzSolver.py:
  this is the main code which calls the other
  classes, it has 4 main functions that produce 
  all the images of the report.

