#!/usr/bin/env python
'''
This program implements Fornberg, 1988
paper for digital differentiators
of arbitrary order.

So, it computes first, second, n derivative along axis 1,2 or 3.
'''


import rsf.api as rsf
import numpy as np
from  scipy.sparse.linalg import cg

functions = {1:goFilter1, 2:goFilter2, 3:goFilter3}
# pars from command line
order  = par.int("order",1) # order of the derivative, default first derivative 
length = par.int("length",5) # filter length, the lengthier the accurate, but also gets costlier 
scale  = par.bool("scale",True) # scales by 1/d^order
axis   = par.int("axis",1) # apply differentiator along axis, default is fast axis

f = Fornberg_filter(length,order)
try:
  functions[axis](Fin,Fout,f,order,scale)
except:
  import sys
  print >> sys.stderr, '========== sfnderiv ERROR ==========='
  print >> sys.stderr, 'Error: valid axis values are 1,2 or 3'
  print >> sys.stderr, '====================================='
  sys.exit(1)

Fin.close()
Fout.close()
