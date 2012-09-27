from jarray import * # Jython module provides function zeros used below
from lab0 import *
from lab0.Dsp import *

print "Jython-Java"
maxtime = 2.0
n = 1001
a = 0.99
x = zeros(n,'f')
y = zeros(n,'f')
x[0] = x[(n-1)/2] = x[n-1] = 1.0
sw = Stopwatch()
sw.start()
nsmooth = 0
while sw.time()<maxtime:
  smooth(a,x,y)
  nsmooth += 1
sw.stop()
print "nsmooth =",nsmooth
print "   mean =",mean(y)
print "   time =",sw.time()
print " mflops =",6.0e-6*n*nsmooth/sw.time()
