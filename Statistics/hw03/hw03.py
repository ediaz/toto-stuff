#!/usr/bin/env python

import numpy as np


Pxy = ([0.15, 0.10, 0.10, 0.10, 0.05], [0.05, 0.08, 0.14, 0.08, 0.05],[0.01, 0.01, 0.02, 0.03, 0.03])


nx = len(Pxy[0])
ny = len(Pxy)


Px = np.zeros(nx)
for i1 in range(nx):
  for i2 in range(ny):
    Px[i1] += Pxy[i2][i1]


Py = np.zeros(ny)
for i2 in range(ny):
  for i1 in range(nx):
    Py[i2] += Pxy[i2][i1]

# 1.a
print 'Px:', Px, 'sum: ',sum(Px)
print 'Py:', Py, 'sum: ',sum(Py)




I = np.zeros((ny,nx))

for i2 in range(ny):
  for i1 in range(nx):
    I[i2][i1]= Px[i1]*Py[i2]


# 1.b
check = (I==Pxy)

# X and Y are not independent
X = range(nx)
Y = range(ny)

def mean(Px):
  mean=0.0
  for i in range(len(Px)):
    mean +=Px[i]*i
  return mean


def var(Px):
  mx = mean(Px)
  var = 0.0  
  for i in range(len(Px)):
    var +=Px[i]*(i-mx)**2
  return var
  

def corr(Pxy,Px,Py):
  nx = len(Px)
  ny = len(Py)

  mx = mean(Px)
  my = mean(Py)

  corr = 0.0 

  for i1 in range(nx):
    for i2 in range(ny):
      corr += (i1-mx)*(i2-my)*Pxy[i2][i1]

  corr = corr/np.sqrt(var(Px)*var(Py))
  return corr









meanx = mean(Px)
varx  = var(Px)

meany = mean(Py)
vary  = var(Py)

# 1.c
print 'mean x:',meanx,'var x:', varx
print 'mean y:',meany,'var y:', vary


Corr_xy  = corr(Pxy,Px,Py)

# 1.d
print 'correlation :', Corr_xy

# 1.e
print 'Conditional probability: P(X | Y=0) = P[X][Y=0]/P[Y=0]'

Px_y0 = np.zeros(nx) 
for i1 in range(nx):
  Px_y0[i1]=Pxy[0][i1]/Py[0]

print Px_y0
print "mean of P(X | Y=0) = ",mean(Px_y0)


#1.f find the mean of D = 5Y +X  => E(D) = 5*mu_y + mu_x 

mu_D = 5*meany +meanx

print "mean D: ",mu_D 




#1.g Find the prob of D<= 7

Pd_7 = 0.0

for i1 in range(nx):
  for i2 in range(ny):
    if(i2*5+i1 <= 7 ):
      Pd_7 += Pxy[i2][i1]

print "prob of D<= 7 :", Pd_7



print 'average number of inspections: ',int(1/Pd_7)





#########
# exercise 2
########

#binomial distribution

n = 600
p = 0.01 

def choose(n,k):
  from scipy.misc import comb
  return(comb(n,k))

def bin_pmf(n,p,k):
  x = choose(n,k) * p**k * (1-p)**(n-k)
  return x


def bin_cdf(n,p,k):
  suma=0.0
  for i in range(k+1):
    suma += bin_pmf(n,p,k)
  return suma

p_bin = 1.0 - bin_cdf(n,p,5)


# Poisson approximation

def poisson_pmf(mu,k):
  from math import exp
  from math import factorial as fact
  x = exp(-mu) * mu**k/fact(k)
  return x

def poisson_cdf(mu,k):
  suma=0.0
  for i in range(k+1):
    suma += poisson_pmf(mu,k)
  return suma


p_poisson = 1.0 - poisson_cdf(n*p,5)



# Exponential approximation
sigma = np.sqrt(n*p*(1-p))
mu    = n*p

from scipy.stats import norm

p_gauss = 1.0 - norm.cdf((5-mu)/sigma)
print 'p binary:',p_bin,'p poisson:',p_poisson,'p_gauss:',p_gauss
