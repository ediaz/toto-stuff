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

n = 1001
maxtime = 2.0 
a =0.99 

x= np.zeros(n, dtype=float)

x[0]=x[int(n/2)]= x[n-1]= 1.0 

nsmooth=0
y= np.zeros(n, dtype=float)

sw.start()
while sw.time() < maxtime :
    dsp.smooth(a,x,y)
    nsmooth+=1
sw.stop()

print "nsmooth = %d "%nsmooth
print "   mean = %9.7f "%dsp.mean(y)
print "   time = %9.7f "%sw.time()
print " mflops = %f "%((6.0e-6*n*nsmooth/sw.time()))



if plot:
    ll = plt.plot(range(n),y)
    plt.show()
