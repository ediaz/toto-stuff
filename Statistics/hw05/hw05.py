#!/usr/bin/env python

import numpy as np
from scipy import stats

from math import *



def bootplanedistance(R):

  mx, sx = 1.0, 0.05
  my, sy = 2.0, 0.1
  dm = sqrt(mx*mx + my*my)
  d = np.zeros(R)
  dnp = d
  Xi = [1.015,0.933,1.034,1.081,0.965,1.043,1.063,0.920,0.928,1.029]
  Yi = [1.981,2.073,1.941,2.218,1.986,2.011,2.107,2.006,1.990,1.917]
  n = len(Xi)
  dm2 = sqrt( np.mean(Xi)**2+np.mean(Yi)**2)
  for i in range(R):
    #parametric
    X =np.random.normal(mx,sx)  
    Y =np.random.normal(my,sy)  
    d[i] = sqrt(X*X +Y*Y)
    # nonparametric
    I = np.floor(n*np.random.uniform(0,1,n))
    xsample = [Xi[int(j)] for j in I]
    ysample = [Yi[int(j)] for j in I]
    XX = np.mean(xsample) 
    YY = np.mean(ysample)
    dnp[i] = sqrt(XX*XX+YY*YY)

  pbias = np.mean(d) -dm
  pstd = np.std(d)

  npbias = np.mean(dnp) -dm2
  npstd = np.std(dnp)

  return [pbias,pstd,npbias,npstd]





print bootplanedistance(5000)

Xi = [1.015,0.933,1.034,1.081,0.965,1.043,1.063,0.920,0.928,1.029]

print np.mean(Xi),np.std(Xi)







def d2x(x,y):
  d2 = (x**2+y**2)**(-0.5) -x*x/(pow(x**2+y**2,3/2))
  return d2




def d2y(x,y):
  d2 = (x**2+y**2)**(-0.5) -y*y/(pow(x**2+y**2,3/2))
  return d2

def bias(x,y,sx,sy):
  bias=0.5*(d2x(x,y)*sx*sx + d2y(x,y)*sy*sy)
  return bias




def dy(x,y):
  dy = y*pow(x*x+y*y,-0.5)
  return dy


def dx(x,y):
  dx = x*pow(x*x+y*y,-0.5)
  return dx


def s2(x,y,sx,sy):
  s2 = pow(dx(x,y),2)*sx*sx + pow(dy(x,y),2)*sy*sy 
  return sqrt(s2)



print "bias:",bias(1.0,2.0,0.05,0.1)

print "std:",s2(1.0,2.0,0.05,0.1)




sigmas=[0.1,0.1,1,1,10,1]
partials= [ -0.249,0.199,-0.00199,0.00199,0.000825,0.000332]
s2=0.0

for i in range(len(sigmas)):
  s2+= pow(partials[i],2)*pow(sigmas[i],2)



print "std2:",sqrt(s2)
