import numpy as np
from scipy import *
from scipy.sparse import *
from scipy.integrate import simps
from scipy.sparse import dia_matrix,linalg
from numpy import linalg as nplin

def integrate(f,a,b,N):
  x = np.linspace(a,b,N)
  nf = len(f)
  fx = f[0]
  y =  x*0
  for i in range(len(x)):
    y[i] = fx(x[i])

  if nf>1:
    for ifn in range(1,len(f)):
      fx = f[ifn] 
      for i in range(len(x)):
        y[i] *= fx(x[i])
  return simps(y, x)





class fem1:
  def __init__(self,a,b,g,sigma,N,DN=False,MIX=False,gamma=1.,upa=1.,upb=1.):
    ''' 
    here, we have N-2 unkowns
    '''
    self.a = a # lower limit of interval
    self.b = b # upper limit of interval
    self.g = g # pointer to function
    self.sigma = sigma # conductivity function
    self.grd = grid1d(a,b,N)
    self.nu = N-2
    self.DN = DN
    self.MIX = MIX
    if DN : self.nu += 1 
    if MIX: 
      self.nu += 2 
      self.gamma = gamma
      self.upa = upa # u prime @ a
      self.upb = upb # u prime @ b

    
    self.xi = self.grd.xi
    

  def apply(self):
    self.getA()
    self.getF()
    self.ui = np.zeros(self.nu)
    return linalg.spsolve(self.A,self.F)

  def getA(self):
    # gets stiffness matrix
    n = self.nu
    grd  = self.grd
    # main diagonal construction
    md = np.zeros(n)
    ld = np.zeros(n)
    ud = np.zeros(n)
    xi = grd.xi
    r = n
    o = 0
    if self.DN: r = n-1
    if self.MIX: 
      r = n - 1
      ih1 = pow(grd.ih[1],2)
      md[0] =  ih1*integrate([self.sigma],xi[0],xi[1],11)
      ld[0] =  -ih1*integrate([self.sigma],xi[0],xi[1],11)
      ud[0] = ld[0]
      o = 1
    
    for i in range(o,r):
      ii = i+1
      if self.MIX:
        ii = i 
      int_i1 = integrate([self.sigma],xi[ii-1],xi[ii  ],11)
      int_i2 = integrate([self.sigma],xi[ii  ],xi[ii+1],11)
      ih1 = pow(grd.ih[ii],2)
      ih2 = pow(grd.ih[ii+1],2)
      md[i]   = int_i1*ih1 + int_i2*ih2 
      if (i<r-1):
        ud[ii] = -int_i2*ih2
      if i>0:
        ld[i-1] = -int_i1*ih1

    if self.DN:
      md[r  ]   = int_i2*ih2
      ld[r-1]   = -int_i2*ih2
      ud[n-1]   = -int_i2*ih2 
  
    if self.MIX:
      int_i1 = integrate([self.sigma],xi[n-2],xi[n-1],11)
      ih1    = pow(grd.ih[n-1],2)
      md[r]  = int_i1*ih1
      md[r] += 1.
      ud[r-1] = -integrate([self.sigma],xi[n-3],xi[n-2],11)*pow(grd.ih[n-2],2)
      ud[r]  = -int_i1*ih1
      ld[r-1]= ud[r]
    self.A = dia_matrix(([ld,md,ud],array([-1,0,1])),shape=(n,n))
    self.A = self.A.tocsr()

  ################################################################
  def getF(self):
    grd  = self.grd
    n = self.nu
    self.F = np.zeros(n)
    r = n
    o = 0
    if self.DN or self.MIX: r = n-1
  
    if self.MIX: 
      o=1
      self.i = 0
      self.F[0] = integrate([self.g,self.phiright],grd.xi[0],grd.xi[1],11)
      self.F[0] -= self.upa*self.sigma(self.a)

    for i in range(o,r):
      ii = i+1
      if self.MIX:
        ii = i

      ll = grd.xi[ii-1]
      ul = grd.xi[ii+1]
      self.i = ii 
      int_i  =  integrate([self.g,self.phileft],ll,grd.xi[ii],11)
      int_i1 =  integrate([self.g,self.phiright],grd.xi[ii],ul,11)
      self.F[i] = int_i+int_i1
    if self.DN or self.MIX: 
      self.F[i+1] = int_i1   
    
    if self.MIX:
      self.F[i+1] += self.sigma(self.b)*self.gamma 



  def phileft(self,x):
    i = self.i 
    grd = self.grd
    a = 1./grd.h[i]
    b = -grd.xi[i-1]/grd.h[i]
    return (a*x+b)


  def phiright(self,x):
    i = self.i 
    grd = self.grd 
    a = -1./grd.h[i+1]
    b = grd.xi[i+1]/grd.h[i+1]
    return (a*x+b) 



##############################################################
class grid1d:
  '''
  this class creates the grid, as it is now
  it only creates a regular mesh. Nevertheless,
  the code is written such that it would accept
  an arbitraty mesh such that:
  x0 < x1 < x2 <...<xn < xn+1
  
  where n is the number of intervals in \Omega =(a,b)
  '''
  def __init__(self,a,b,N,Constant=True):
    self.a = a
    self.b = b
    self.n = N
    if Constant:
      self.xi = np.linspace(self.a,self.b,self.n)
    self.hi()

  def hi(self):
    self.h = np.zeros(self.n)
    self.ih = np.zeros(self.n)
    self.h[0] = 0.
    for i in range(1,self.n):
      self.h[i] = self.xi[i]-self.xi[i-1]
      self.ih[i] = 1./self.h[i]

  def geti(self,x):
    j = -1
    for i in range(self.n):
      if x < self.xi[i]:
        j = i
        break
    return j 

################################################################
class ux:
  '''
  given a FEM ui(xi) solution, get the interpolated
  continous function u(x)
  '''
  def __init__(self,ui,grd,N,f):
    self.ui  = ui
    self.grd = grd
    self.f = f
    xi = grd.xi
    a = xi[0]
    b = xi[len(xi)-1]
    self.x = np.linspace(a,b,N)

  def L2(self):
    self.eval()
    diff = (self.f(self.x) - self.u)
    for i in range(len(diff)):
      diff[i] *= diff[i]
    l2 = simps(diff,self.x)
    return np.sqrt(l2)

  def eval(self):
    n = len(self.x)
    grd = self.grd
    u = np.zeros(n)
    for i in range(n):
      x = self.x[i]
      j = grd.geti(x)
      wi = self.phileft(x,j)
      u[i] = self.ui[j]*wi +(1.-wi)*self.ui[j-1] 
    self.u = u


  def phileft(self,x,i):
    grd = self.grd
    a = 1./grd.h[i]
    b = -grd.xi[i-1]/grd.h[i]
    return (a*x+b)

  def phiright(self,x,i):
    grd = self.grd 
    a = -1./grd.h[i+1]
    b = grd.xi[i+1]/grd.h[i+1]
    return (a*x+b) 


