#!/usr/bin/env python

import numpy as np
import numpy as np
from matplotlib import pyplot as plt
import qqplot as qq
import scipy.stats as stats

import math

sulphur=np.array([52.7,53.9,41.7,71.5,47.6,55.1,
                  62.2,56.5,33.4,61.8,54.3,50.0,
                  45.3,63.4,53.9,65.5,66.6,70.0,
                  52.4,38.6,46.1,44.4,60.7,56.4])


qq.qqplots(sulphur,title='',xlabel='z-score',ylabel='normal vaue')
plt.show()


print "==============================="
print "          Problem 4            "
shap = stats.shapiro(sulphur)
print "shapiro output\n",shap
print "mean shulpur = ",np.mean(sulphur)
print "std shulpur = ",np.std(sulphur)

print np.mean(sulphur),"+/-",stats.norm.ppf(1-0.05/2)*np.std(sulphur)/math.sqrt(24)


print "Z=",(54.3333-50)/(np.std(sulphur)/math.sqrt(24))

za= abs(stats.norm.ppf(0.05/2))

print "Z_(alpha=0.05) = ",za
print "P(Z>za) =",stats.norm.sf(za)


########## Problem 5 ##############

print "==============================="
print "          Problem 5            "
Ip1 = 3./171. 
Ip2 = 8/55. 

Ep1 = 0/123.
Ep2 = 12/122.


def ciprop(p1,p2,n1,n2,alpha=0.05):
  mu = p1-p2
  za = abs(stats.norm.ppf(alpha/2))
  sq = math.sqrt((p1*(1-p1)/n1)+(p2*(1-p2))/n2)
  ci = [mu -za*sq,mu +za*sq]

  return ci 
  


#consistency amog studies P1:
ci1 = ciprop(Ip1,Ep1,171,123)

print "ci for infections using condoms:",ci1
ci2 = ciprop(Ip2,Ep2,55,122)

print "ci for infections no condoms:",ci2

ci3 = ciprop(Ip1,Ip2,171,55)
print "ci testing H0 = no correlation between using or not condoms (I):",ci3

ci4 = ciprop(Ep1,Ep2,123,122)
print "ci testing H0 = no correlation between using or not condoms (E):",ci4

########################################################################




print "==============================="
print "          Problem 2            "


def dFdto(x,to,NMO):
  dx = -0.5*x*math.sqrt(1/(2*NMO))*pow(to,-1.5)
  return dx

def dFdNMO(x,to,NMO):
  dx = dFdto(x,NMO,to)
  return dx


def mu(x,to,NMO):
  mu = math.sqrt(x*x/(2*NMO*to))
  return mu


def sigma(x,to,NMO,sto,sNMO):
  d1 = dFdto(x,to,NMO)
  s1 = sto

  d2 = dFdNMO(x,to,NMO)
  s2 = sNMO

  s3 = d1*d1*s1*s1 + d2*d2*s2*s2

  return math.sqrt(s3)


print "mu = ",mu(200.,1.,0.005)
print "sigma = ",sigma(200.,1.,0.005,0.005,0.001)

print " DF2 = %f"%dFdNMO(200.,1.,0.005)


print "==============================================="
print "             Problem 6                         "
n1 = 1997
n2 = 906
n3 = 904
n4 = 32

a = n1 + n2 +n3+n4
b = -(n1-2*n2-2*n3-n4)
c = -2*n4

x = (-b +math.sqrt(b*b -4*a*c))/(2*a)
print "roots of theta: ",x





p1 = 0.25*(2+x)*a
p2 = 0.25*(1-x)*a
p3 = 0.25*(1-x)*a
p4 = 0.25*x*a


obs = [n1,n2,n3,n4]
pre = [p1,p2,p3,p4]
xx = 0.0 
for i in range(len(obs)):
  print i," &", obs[i]," &", "%5.2f"%pre[i], "\\\\"
  num = (obs[i]-pre[i])
  num *= num
  xx += num/pre[i]

print xx

p =  stats.chi2.sf(xx,2)
print p




print "==============================================="
print "             Problem 7                         "

h1 = [0,1,9,10]

alpha =0 
for x in (h1):
  alpha += stats.binom.pmf(x,10,0.6)

print alpha
