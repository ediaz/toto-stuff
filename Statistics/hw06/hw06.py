#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
import qqplot as qq
import scipy.stats as stats

import math






#problem 1:
s1,n1 = 0.7,100
s2,n2 = 0.6,120
s3,n3 = 0.8,130

ss = 0.25*(s1*s1+s2*s2)+s3*s3
ss = ss/(n1+n2+n3-3)
print ss, math.sqrt(ss)


tmp = stats.norm.ppf(1-0.05/2)*math.sqrt(0.25*(s1*s1/n1 +s2*s2/n2)+s3*s3/n3)
print (-1.2-tmp,-1.2+tmp)


print stats.norm.ppf(1-0.05/2)

#problem 2

whitecollar=np.array([28.57,20.10,69.05,65.40,29.59,44.82,77.37,24.67,65.01,9.99,12.20,22.55,14.30,31.79,11.6,68.47,42.64,16.70,86.27,76.73])

#qq.qqplots(whitecollar,title='',xlabel='z-score',ylabel='normal vaue')
#plt.show()


sd = np.std(whitecollar)
xbar =  np.mean(whitecollar)
n= len(whitecollar)
tmp = abs(stats.t.ppf(0.05/2,n-1) * sd/math.sqrt(n))

ci = (xbar -tmp,xbar+tmp)

print ci


#problem 3
ground = np.array([4.446,3.99,3.73,3.29,4.82,6.71,4.61,3.87,3.17,4.42,3.76,3.3]) 
sateli = np.array([4.08,3.94,5.00,5.2,3.92,6.21,5.95,3.07,4.76,3.25,4.89,4.80])

dif = [ground[i]-sateli[i] for i in range(len(ground))]

mdif = np.mean(dif)
mstd = np.std(dif)
n = len(ground)
tmp  =  abs(stats.t.ppf(0.05/2,n-1) * mstd/math.sqrt(n))

ci = (mdif - mstd,mdif+mstd)
print ci


#problem 7.3
m = 6.3
sd = 4.57
tmp = stats.norm.ppf(1-0.05/2)*sd/math.sqrt(96)
print (m -tmp,m+tmp)

