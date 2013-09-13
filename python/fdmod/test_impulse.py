import numpy as np
from rsflikeLibrary import *
import plot as plt



# from Francesco's modeling code:
C1 = +0.598144
C2 = -0.039876
C3 = +0.004785
C4 = -0.000348 


def DerivForward(inp):
  n = len(inp)
  out = np.zeros(n)
  for iz in range(4,n-4,1):
    out[iz] = 2*C4*(inp[iz+4]-inp[iz-3])+\
              2*C3*(inp[iz+3]-inp[iz-2])+\
              2*C2*(inp[iz+2]-inp[iz-1])+\
              2*C1*(inp[iz+1]-inp[iz  ])
  out[0] = out[1] = out[2] = out[3] = out[4]
  out[n-1] = out[n-2] = out[n-3] = out[n-1]
  return out



def DerivBackward(inp):
  n = len(inp)
  out = np.zeros(n)
  for iz in range(4,n-4,1):
    out[iz] = 2*C4*(inp[iz+3]-inp[iz-4])+\
              2*C3*(inp[iz+2]-inp[iz-3])+\
              2*C2*(inp[iz+1]-inp[iz-2])+\
              2*C1*(inp[iz  ]-inp[iz-1])
  out[0] = out[1] = out[2] = out[3]
  out[n-1] = out[n-2] = out[n-3] = out[n-4]=out[n-5]
  return out




x = np.linspace(0,10.,101)
inp = x**2.
forw = DerivForward(inp)
bacw = DerivBackward(forw*(10))

rsfforw = rsf1darray(0.0,0.01,101)
rsfforw.put(forw*10.)

rsfbacw = rsf1darray(0.0,0.01,101)
rsfbacw.put(bacw*(10.))
#plt.plot1d(rsfbacw)
#plt.plot1d(rsfforw)
#
#rsffft = rsffft1(rsfforw)
#rsfforwfft = rsffft.fft()
#rsfforwfft.rsf = np.abs(rsfforwfft.rsf)
#plt.plot1d(rsfforwfft)
#rsfforwfft.rsf  = (rsfforwfft.grid()*2.*np.pi)**2
#plt.plot1d(rsfforwfft)
#
#inp = np.zeros(101)
#inp[50] = 1.0
#
#forw = DerivForward(inp)
#
#for i in range(45,56):
#  print "%13.8f%4d"%(forw[i],i)
#
#
import difference as diff

inp *= 0.

inp[50] = 1.



d = diff.DifferenceFilter(2,9)

f = d.filter()
print f
inp2 = np.convolve(inp,f,'same')
rsf = rsf1darray(0.0,0.1,101)
rsf.rsf = inp2
rsfft = rsffft1(rsf)
fft = rsfft.fft()
fft.rsf = np.abs(fft.rsf)
plt.plot1d(fft,fn=1,fname=None)

inp *= 0.

inp[50] = 1.
inp2 = DerivBackward(DerivForward(inp))
rsf = rsf1darray(0.0,0.1,101)
rsf.rsf = inp2
rsfft = rsffft1(rsf)
fft = rsfft.fft()
fft.rsf = np.abs(fft.rsf)
plt.plot1d(fft,fn=1,fname=None)





plt.show()
#
#
#
#print x,a
