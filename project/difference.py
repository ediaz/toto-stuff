import numpy as np
import matplotlib.pyplot as plt


class DifferenceFilter:
  def __init__(self,order=2,n=11,h=1.):
    self.order  = order # derivative order, default is second order
    self.n = n # lenght of the filter (is centered)
    self.h = h # interval of discretization for scaling the filter
    self.ih = (1./h)**order
  
  def filter(self):
    c = self.Fornberg_filter()
    self.ir = c[self.order]
    return self.ir*self.ih
  ######################## Fornberg ################################
  def Fornberg_filter(self):
    x = np.arange(-int(self.n/2),self.n/2+1,1)
    z = 0.0
    k = self.order
    n = self.n 
    m = n -1
    c1 = 1
    c4 = x[0] - z
    c = np.zeros((n, m+1))
    c[0,0] = 1
    for i in xrange(1, n):
      mn = min(i, m)
      c2 = 1
      c5 = c4
      c4 = x[i] - z
      for j in xrange(i):
        c3 = x[i] - x[j]
        c2 = c2*c3
        if j == i-1:
          for v in xrange(mn, 0, -1):
            c[i,v] = c1*(v*c[i-1,v-1] - c5*c[i-1,v])/c2
          c[i,0] = -c1*c5*c[i-1,0]/c2

        for v in xrange(mn, 0, -1):
          c[j,v] = (c4*c[j,v] - v*c[j,v-1])/c3
        c[j,0] = c4*c[j,0]/c3
      c1 = c2
    return c.T[0:k+1]  

  ########################   Graphic QC    ##########################
  def showfilter(self,fn=1):
    a = np.zeros(200)
    a[100] = 1.
    n = self.ir
    c = np.convolve(a,self.ir,mode='same')
    fft_c = np.fft.fftshift(np.fft.fft(c))
    y = np.abs(fft_c)
    x = np.arange(-0.5,0.5,1./len(fft_c))
    yy = np.abs((x*2.*np.pi)**self.order)
    fig = plt.figure(fn)
    ax = fig.add_subplot(111)
    p1 = ax.plot(x,y)
    p2 = ax.plot(x,yy)
    ax.set_xlim(-0.5,0.5) 
    ax.set_xlabel('f(cycles/sample)')

    ax.legend(('%s method impulse response N=%d, order=%d'%(self.method,len(self.ir),self.order),\
'theoretical response: $\left|(2i\pi f)^%d\\right|$'%self.order),'upper center')




