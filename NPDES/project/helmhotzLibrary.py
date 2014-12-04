import numpy as np
from difference import DifferenceFilter as Diff
from rsflikeLibrary import *
import scipy.sparse as sparse
from scipy.sparse import linalg as splinalg
from copy import copy
import matplotlib.pyplot as plt

class helmhotz:
  def __init__(self,rsfvel,rsfwav,sx,sz,nw=None,\
               free=False,\
               N=11,abc=True,verb=False,f=10.):
    self.slo = copy(rsfvel)
    self.slo.prop = 'slo^2'
    self.slo.Punit = 's^2/km^2'
    self.slo.rsf = 1./(self.slo.rsf)**2
    self.n = self.slo.n2*self.slo.n1

    if rsfwav != None:
      wobj = rsffft1(rsfwav)
      self.rsfwav = wobj.fft()
      self.sx,self.sz = sx,sz
      if verb:
        print 'frequency grid:'
        self.rsfwav.geom()
    if verb:
      print 'spatial grid:'
      self.slo.geom()

    if nw == None: 
      self.nw = self.rsfwav.n1
    else:
      self.nw = nw

    self.wavefield = np.zeros((self.nw,self.slo.n2,self.slo.n1),dtype=complex)
    self.N = N # number of samples for derivative filters
    fd2 = Diff(order=2,n=N)
    self.fd2 = fd2.filter()
    fd1 = Diff(order=1,n=N)
    self.fd1 = fd1.filter()
    self.pml = abc
    self.free = free
    self.f = f 

  def solve(self):
    '''
    Here, I solve the linear system of equations
    '''
    isx = int((self.sx - self.slo.o2)/self.slo.d2)
    isz = int((self.sz - self.slo.o1)/self.slo.d1)
    isn = self.slo.n1*isx+isz
    for iw in range(self.nw):
      b = np.zeros(self.n,dtype=complex)
      f = self.rsfwav.o1 +self.rsfwav.d1*iw
      iw2 = iw
      if self.nw == 1: 
        f = self.f
        iw2 = int((f -self.rsfwav.o1)/self.rsfwav.d1)
      print 'iw =',iw2,'f=',f,self.rsfwav.rsf[iw2]
      w = 2.*np.pi*f
      
      b[isn] -= self.rsfwav.rsf[iw2]
      A = self.lhs(w)
      u = splinalg.spsolve(A,b)
      self.wavefield[iw] = np.reshape(u,(self.slo.n2,self.slo.n1)) 
    self.A = A
    return self.wavefield

  def returnA(self):
    return self.A

  def solveE(self,w,b):
    A = self.lhs(w)
    u = splinalg.spsolve(A,b)
    u = np.reshape(u,(self.slo.n2,self.slo.n1))
    return u 

  def lhs(self,w):
    self.buildpml(w)
    return self.k2(w)+self.laplace() +self.grads()   

  def buildpml(self,w):
    ''' 
    This routine  sets the PML class object
    '''
    abc = self.pml
    if w== 0.0 : abc= False
    if abc:
      pmlx = pml(self.slo.o2,self.slo.d2,self.slo.n2)
      pmlz = pml(self.slo.o1,self.slo.d1,self.slo.n1,free=self.free)
      self.px = pmlx.wpml(w)
      self.dpx = pmlx.dwpml(w)
      self.pz = pmlz.wpml(w)
      self.dpz = pmlz.dwpml(w)
      pmlx = pmlz = None
    else:
      self.px = np.ones(self.slo.n2)
      self.dpx = np.zeros(self.slo.n2)
      self.pz = np.ones(self.slo.n1)
      self.dpz = np.zeros(self.slo.n1)

  def k2(self,w): 
    '''
    This does the wavenumber scaling matrix
    '''     
    k2 = np.reshape(w*w*self.slo.rsf,(self.n))
    return sparse.dia_matrix((k2,[0]),shape=(self.n,self.n))

  def laplace(self):
    ''' 
    This function returns the laplacian matrix operator:
    rx^2*d_xx + rz^2*d_zz
    '''
    f = self.fd2
    nx = self.slo.n2
    nz = self.slo.n1
    def Bh(n,r):
      nonzero = np.ones((self.N,n))
      offsets = np.array(range(self.N/2,-self.N/2,-1))
      for i in range(self.N):
        nonzero[i] *= self.fd2[i]
      return sparse.spdiags(r*r,0,len(r),len(r))*sparse.spdiags(nonzero,offsets,n,n)
    #
    def I(n):
      return sparse.spdiags(np.ones(n),[0],n,n)

    Bx = (1./self.slo.d2**2)*Bh(nx,self.px); Ix = I(nx)
    Bz = (1./self.slo.d1**2)*Bh(nz,self.pz); Iz = I(nz)
    return sparse.kron(Bx,Iz) + sparse.kron(Ix,Bz)


  def grads(self):
    ''' 
    This function returns the sum of the derivatives operator:
    rx*r_x*d_xx + rz*r_z*d_z
    '''    
    nx = self.slo.n2
    nz = self.slo.n1
    def Bh(n,r):
      nonzero = np.ones((self.N,n))
      offsets = np.array(range(self.N/2,-self.N/2,-1))
      for i in range(self.N):
        nonzero[i] *= self.fd1[i]
      return sparse.spdiags(r,0,len(r),len(r))*sparse.spdiags(nonzero,offsets,n,n)
    def I(n):
      return sparse.spdiags(np.ones(n),[0],n,n)
    Bx = (1./self.slo.d2)*Bh(nx,self.px*self.dpx); Ix = I(nx)
    Bz = (1./self.slo.d1)*Bh(nz,self.pz*self.dpz); Iz = I(nz)
    return sparse.kron(Bx,Iz) + sparse.kron(Ix,Bz)
  

################################################################################
class pml:
  def __init__(self,o,d,n,c=1.0,free=False,power=2,so=40.):
    self.o = o
    self.d = d
    self.n = n
    self.free = free
    self.power = power 
    self.so = so
    self.c = c
    self.xmax = o+(n-1)*d
    if free:
      self.weight = self.wfree
    else:
      self.weight = self.wnfree

  def wfree(self,x):
    return np.abs(x-self.xmax)  

  def wnfree(self,x):
    d1 = (x - self.o)
    d2 = -(x - self.xmax)
    return np.abs(min(d1,d2))

  def wpml(self,w):
    self.lm = self.c*2.*np.pi/(w*self.d)  # width of the pml = 1 wavelenght (in number of samples)
    self.w = w
    x = np.linspace(self.o,self.xmax,self.n)
    wght = np.zeros(self.n,dtype=complex)
    id1 = 1./self.d
    for i in range(self.n):
      xw =  self.weight(x[i])
      ixw = xw*id1
      if ixw <= self.lm:
        d = (self.lm - ixw)
        s = self.so*(d)**self.power/(self.lm**self.power) 
      else:
        s = 0.  # remains equal to 0
      wght[i] = 1./(1+1.0j*s/w)
    return wght

  def dwpml(self,w):
    self.lm = self.c/w  # width of the pml = 1 wavelenght
    self.w = w
    x = np.linspace(self.o,self.xmax,self.n)
    wght = np.zeros(self.n,dtype=complex)
    if self.power%2. == 0:
      p = self.power +1
    else:
      p = self.power
    for i in range(self.n):
      xw =  self.weight(x[i])
      if xw <= self.lm:
        d = (self.lm-xw)
        s = self.so*d**self.power 
        sd = self.so*self.power*(d**(self.power-1))
        wi = 1./(1.+1.0j*s/w)
      else:
        wi = 0.
        sd = 0.0
      wght[i] = wi*wi*(1.0j/w)*sd*np.sign(x[i]-0.5*(self.o +(self.n-1)*self.d))**(p)
    return wght



