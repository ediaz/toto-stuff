import rsf.api as rsf
import numpy as np
import sys



def get_axis(File,axis):
  o = File.float("o%d"%axis)
  d = File.float("d%d"%axis)
  n = File.int("n%d"%axis)
  return o,d,n

def put_axis(File,axis,o,d,n):
  File.put("o%d"%axis,o)
  File.put("d%d"%axis,d)
  File.put("n%d"%axis,n)  

par = rsf.Par()

Fin = rsf.Input()
Fout = rsf.Output()
Fp = rsf.Input("pen")

o1,d1,n1 = get_axis(Fp,1)

print n1
penalty = np.zeros((n1),'f')
Fp.read(penalty)


o1,d1,n1 = get_axis(Fin,1)
o2,d2,n2 = get_axis(Fin,2)
o3,d3,n3 = get_axis(Fin,3)

penaltyExt = np.zeros((n2,n1),'f')
panel  = np.zeros((n2,n1),'f')

for i2 in range(n2):
  for i1 in range(n1):  
    penaltyExt[i2,i1] = penalty[i2]


for i3 in range(n3):
  Fin.read(panel)
  Fout.write(panel*penaltyExt)


Fin.close()
Fout.close()
Fp.close()

