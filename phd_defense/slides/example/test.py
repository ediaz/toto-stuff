import rsf.api as rsf
import numpy as np
import sys


def corr(f,g,c,lags,adj=False):
  nlags = len(lags)
  for ci in range(nlags):
    l = lags[ci] 
    for i in range(len(f)):
      if i+l < len(f) and i+l >= 0 and i-l >=0 and i-l <len(f):
        if not adj:
          c[ci ] += f[i-l]*g[i+l]
        else:
          f[i-l] += c[ci ]*g[i+l]


def put_axis(File,axis,o,d,n):
  File.put("o%d"%axis,o)
  File.put("d%d"%axis,d)
  File.put("n%d"%axis,n)  

par = rsf.Par()

Fin = rsf.Input()
Fout = rsf.Output()
Fg = rsf.Input("g")


nl = par.int("nlag")
 
lags = range(-nl,nl+1,1)


adj = par.bool("adj",False)

n1 = Fg.int("n1")
g = np.zeros(n1,'f')
Fg.read(g)

f = np.zeros(n1,'f')
if adj:
  nl = int(Fin.int("n1")/2)
else:
  Fin.read(f)
  
c = np.zeros(nl*2+1,'f')

if adj:
  Fin.read(c)

corr(f,g,c,lags,adj)


if adj:
  put_axis(Fout,1,0,Fg.float("d1"),n1)
  Fout.write(f)
else:
  put_axis(Fout,1,-nl*2*Fg.float("d1"),2*Fg.float("d1"),nl*2+1)
  Fout.write(c)



Fin.close()
Fout.close()
Fg.close()

