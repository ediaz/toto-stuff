import numpy as np
import scipy as sp
from femlib import *
from util   import *
from scipy.integrate import quadrature as quad
import matplotlib.pyplot as plt
import sys


plt.rcParams.update({'font.size': 15,'legend.fontsize': 15})


'''
  PLOTEA EL ERROR UH NO UJ!!!!!!!!!!!!!!!

'''



NI = 1001
pi = np.pi

goPlot = False 

if goPlot:
  Ks = [0]
else:
  Ks = [0,1,2,3,4]



def goP1():
  # problem 1:
  a = 0.0
  b = pi
  L2 = 1. 

  print "Problem 1: Dirichlet BC"
  print "k&   N &      h    &     L2    &     rate"
  for k in Ks:
    N = pow(2,k)*100 +1
    obj = fem1(a,b,g1,s1,N) # sets up the FEM class object
    u1 = obj.apply()             # gets the solution at the free nodes
    u = np.zeros(N)              
    u[1:N-1] = u1                # includes the value at constrained nodes
    interp = ux(u,obj.grd,NI,np.sin) # gets u(x) = \sum_i u_i phi_i(x)
    r = np.log(L2/interp.L2())/np.log(2.)
    print "%d& %4d& %10.8f& %10.8g& %10.8g"%(k,N,(b-a)/(N-1),interp.L2(),r)
    L2 =  interp.L2()
  
  #printAF(obj.A.todense(),obj.F)
  grd = obj.grd
  f = np.sin 
  if goPlot: plot(1,f,interp) 

def goP2():
  #problem 2:
  a = 0.0
  b = pi/4.
  L2 = 1. 
  print "Problem 2: Dirichlet-Neumman BC"
  print "k&   N &      h    &     L2    &     rate"
  for k in Ks:
    N = pow(2,k)*100 +1 
  
    obj = fem1(a,b,g2,s2,N,DN=True)
    u1  = obj.apply()
    u = np.zeros(N)
    u[1:N] = u1
    interp = ux(u,obj.grd,NI,f2)
    r = np.log(L2/interp.L2())/np.log(2.)
    print "%d& %4d& %10.8f& %10.8g& %10.8g"%(k,N,(b-a)/(N-1),interp.L2(),r)

    L2 =  interp.L2()

  grd = obj.grd
  f = f2
  if goPlot: plot(2,f,interp) 



def goP3():
  a = 0.0
  b = 1.0
  L2=1.
  print "Problem 3: Mixed BC"
  print "k&   N &      h    &     L2    &     rate"
  for k in Ks: 
    N = pow(2,k)*100 +1
    obj = fem1(a,b,g3,s3,N,MIX=True,upa=3.,gamma=4.*np.exp(3.))
    u  = obj.apply()
    interp = ux(u,obj.grd,NI,f3)     
    r = np.log(L2/interp.L2())/np.log(2.)
    print "%d& %4d& %10.8f& %10.8g& %10.8g"%(k,N,(b-a)/(N-1),interp.L2(),r)

    L2 =  interp.L2()
    
  grd = obj.grd
  f = f3
  if goPlot: plot(3,f,interp)
  





goP1()
goP2()
goP3()

plt.show()
  
