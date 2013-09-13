#!/usr/bin/env python

import stopwatch  as sw
import dsp as dsp 
import numpy as np 

try: 
    import matplotlib.pyplot as plt
    plot=True
except:
    print "sorry, I won't be able to plot"
    plot=False



print "Python"

n = 2500
a = 0.99

x = np.genfromtxt("SpringDashpot.txt", dtype=None)
#x[0]=x[int(n/2)]= x[n-1]= 1.0 

for i in range(len(x)):
  x[i] += 1.


nsmooth=0
y= np.zeros(n, dtype='Float32')

for i in range(12):
  dsp.smooth2(a,x,x)




if plot:
    ll = plt.plot(range(n),x)
    plt.show()
