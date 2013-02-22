The software in this report was written in Python.
For readability, I splitted the code in three modules, and a main script:

functions_hw01.py: 
  contains the force and boundary functions

plotting_hw01.py:
  contains plotting functions

heatEquation.py:
  contains the C-N class to solve the heat equation

hw01.py:
  main script, loops over sampling h,dt and constructs
  the plots

Dependencies:
  The results of this homework were built using python2.7
  Some external python packadges are required:
    Scipy
    Numpy
    Matplotlib


========================================================================
The report is in the folder report, I included the .tex files
so the report can be reproduced.

To run the program:
  in the command line:
  $ python hw01.py

To reproduce the report:
  $ cd report; pdflatex report.tex; pdflatex report.tex


This instructions should work on an unix based OS (linux, OSX,etc).
In Windows it should still work, but I didn't test it.

