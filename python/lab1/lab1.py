#!/usr/bin/env python

import stopwatch2 as sw
import numpy as np 
import lab1

print "Python-Fortran"

n = 1001
maxtime = 2.0 
a =float(0.99) 

x= np.zeros(n, dtype='float32')

x[0]=x[int(n/2)]= x[n-1]= 1.0 

nsmooth=0
y= np.zeros(n, dtype='float32')

# initialize stopwatch object:
t= sw.stopwatch()


t.start()
while t.time() < maxtime :
    lab1.dsp.smooth_dsp(a,x,y)
    nsmooth+=1
t.stop()

print "nsmooth = %d "%nsmooth
print "   mean = %9.7f "%lab1.dsp.mean_dsp(y)
print "   time = %9.7f "%t.time()
print " mflops = %f "%((6.0e-6*n*nsmooth/t.time()))

