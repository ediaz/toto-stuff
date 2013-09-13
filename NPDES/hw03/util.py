import numpy as np
import scipy as sp
from scipy.integrate import quadrature as quad
import matplotlib.pyplot as plt
import sys

def printAF(A,F):
  n = len(A)
  print 'A:'
  form = '%5.4f'
  for i in range(n):
    for j in range(n):
      if A[i,j] == 0:
        f = '%9d'
      else:
        f = '%9.4f'
      if j == n-1:
        f += '\n'
      sys.stdout.write(f % A[i,j])
  
  print 'F:'
  for i in range(n):
    print F[i]
    
def f1(x):
  return np.sin(x)

def g1(x):
  return (1+x)*np.sin(x) - np.cos(x)

def s1(x):
  return (x+1.)



def f2(x):
  return np.sin(2*x)

def g2(x):
  return 4*np.sin(2*x)

def s2(x):
  if type(x) == np.float64 or type(x) == type(1.):
    f = 1.
  else:
    f = np.ones(len(x))
  return f

def f3(x):
  return np.exp(3.*x)

def g3(x):
  return -9.*np.exp(3.*x)  

def s3(x):
  return s2(x)



def plot(fn,f,interp):
  fig = plt.figure(fn,figsize=(8,8))
  ax1 = fig.add_subplot(211)
  ax1.plot(interp.x,interp.u)  
  ax1.plot(interp.x,f(interp.x))  
  ax1.set_xlabel("x")
  ax1.set_ylabel("u")

  ax2 = fig.add_subplot(212)
  ax2.plot(interp.x,f(interp.x)-interp.u)  
  ax2.set_xlabel("x")
  ax2.set_ylabel("e")
  ax2.ticklabel_format(style='sci',scilimits=(1,10))
  fig.savefig('report/fig/fig%d'%fn)
