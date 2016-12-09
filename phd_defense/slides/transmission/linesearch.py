 #!/usr/bin/env python
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()
Fin = rsf.Input()
Fout = rsf.Output()

Fout.put("n1",1)
Fout.put("o1",0.)
Fout.put('d1',1.)

d1 = Fin.float("d1")
n1 = Fin.int("n1")
o1 = Fin.float("o1")
objt = np.zeros(n1,'f')
alpha = np.zeros(n1,'f')

alpha = np.arange(o1,d1*n1+o1,d1)
Fin.read(objt)

C=np.polyfit(alpha,objt,2)

vertex = -.5*C[1]/C[0]
if C[0] < 0 and objt[2]<objt[0]:
  vertex = alpha[2]
if C[0] < 0 and objt[2] > objt[0]:
  vertex = alpha[1]*0.5

if vertex < 2*alpha[2]: vertex=alpha[2]

print >>sys.stderr, vertex


out = np.zeros(1,'f')
out = np.asarray(vertex,'f')
Fout.write(out)

Fout.close()
Fin.close()
