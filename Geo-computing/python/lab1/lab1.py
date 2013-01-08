#!/usr/bin/env python

import stopwatch2 as sw
import numpy as np 
import lab1

try: 
    import matplotlib.pyplot as plt
    plot=True
except:
    print "sorry, I won't be able to plot"
    plot=False



print "Python-Fortran"
maxtime= 2
n = 2500
a = float(0.9998)

x = np.genfromtxt("SpringDashpot.txt", dtype='float32')
for i in range(len(x)):
  x[i] = float(1.)
nsmooth=0
y= np.zeros(n, dtype='float32')

# initialize stopwatch object:
t= sw.stopwatch()


t.start()
while t.time() < maxtime :
    lab1.dsp.smooth_dsp2(a,x,y)
    nsmooth+=1
t.stop()

print "nsmooth = %d "%nsmooth
print "   mean = %9.7f "%lab1.dsp.mean_dsp(y)
print "   time = %9.7f "%t.time()
print " mflops = %f "%((6.0e-6*n*nsmooth/t.time()))

if plot:
    ll = plt.plot(range(n),y)
    plt.show()
