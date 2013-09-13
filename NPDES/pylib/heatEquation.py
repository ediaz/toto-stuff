import numpy as np
from scipy import *
from scipy.sparse import dia_matrix,linalg

'''
This module provides two solution approaches for the 
1d heat-diffusion equation.

explicit Euler class: 
  solves the 1d heat equation explicitly by finite differences


implicit Euler class: 
  solves the 1d heat equation implicitly using a linear system
  of equations.

  [I + a*dt*Bh] U^m = [U^{m-1} +dt*F^m] 

  A = [I + a*dt*Bh]
  b = [U^{m-1} +dt*F^m] 
  x = U^m

  so 
  A x = b

  A = A^t is SPD
'''


class explicitEuler:
  '''
  solves by FD the heat equation:

  U^{m+1}_{n}-U^{m}_{n}     U^{m}_{n+1} -2U^{m}_{n} +U^{m}_{n-1} 
  --------------------- -a* ------------------------------------  = f^m_n 
          dt                            h^2
  '''
  def __init__(self,f,v,\
               nt=11,ot=0.,dt=1.0,\
               nx=11,ox=0.,dx=1.,a=1.,verb=True):
    ''' 
    init object
    '''
    self.nt = nt 
    self.ot = ot
    self.dt = dt

    self.nx = nx
    self.ox = ox
    self.dx = dx
    self.a  = a
    self.u = np.zeros((nt,nx))
    self.f = f(nt,nx) # pointer to a function with parameters nt,nx
                      # it have to return a nx*nt matrix

    self.v = v(nx)    # pointer to a initial boundary function with 
                      # parameter nx               
                      # it has to return a nx vector 

    self.verb = verb
    self.alpha = self.a*self.dt/(self.dx*self.dx)
    if(self.alpha>0.5):
      print "Warning! alpha=%g>0.5, propagation will be unstable"%self.alpha
      print "dt = %g should be fine"%(0.5*(self.dx*self.dx)/self.a)
      print "change M to %d"%(1+int(1/(0.5*(self.dx*self.dx)/self.a)))

    if not(len(self.f)==nt and len(self.f[0])==nx):
      raise "Error: function returns wrong \
             array shape, I expected a nt by nx matrix"

    if not(len(self.v) == nx):
      raise "Initial boundary vector has wrong dimensions \
             , I expected a vector with nx length"
  
  def solve(self):
    self.u[0,:] = self.v
    p10=int(self.nt/10.)

    for it in range(0,self.nt-1,1):
      for ix in range(1,self.nx-1,1):
        lx = (self.u[it,ix+1]-2.*self.u[it,ix]+self.u[it,ix-1])
        self.u[it+1,ix] = self.u[it,ix] + self.alpha*lx +\
                          self.dt*self.f[it,ix]
      if(it%p10==0 and self.verb): print "%3d%% done"%(it/p10*10)
    return self.u 





class implicitEuler:
  '''
  solves the heat equation implicitly using the folowing linear system:
  [I + a*dt*Bh] U^m = [U^{m-1} +dt*F^m] 
       A x = b
  '''
  def __init__(self,f,v,nt=101,ot=0.,dt=1.0,nx=201,\
                               ox=0.,dx=1.,a=1.,verb=True):
    ''' 
    init object
    '''
    self.nt = nt 
    self.ot = ot
    self.dt = dt

    self.nx = nx
    self.ox = ox
    self.dx = dx
    self.a  = a
    self.u = np.zeros((nt,nx))
    self.f = f(nt,nx) # pointer to a function with parameters nt,nx
                      # it have to return a nx*nt matrix

    self.v = v(nx)    # pointer to a initial boundary function with 
                      # parameter nx               
                      # it has to return a nx vector 
    self.verb = verb

    if not(len(self.f)==nt and len(self.f[0])==nx):
      raise "Error: function returns wrong \
             array shape, I expected a nt by nx matrix"

    if not(len(self.v) == nx):
      raise "Initial boundary vector has wrong dimensions \
             , I expected a vector with nx length"
 

  def solve(self):
    self.u[0,:] = self.v

    A = self._bh()
    p10=int(self.nt/10.)
    
    for it in range(1,self.nt,1):
      b = self._rhs(it)       
      x = linalg.spsolve(A,b) # linear solver for sparse A
      for ix in range(1,self.nx-1,1):
        self.u[it,ix] = x[ix-1]
      if(it%p10==0 and self.verb): print "%3d%% done"%(it/p10*10)
    return self.u

  ######## Private functions within the object #########
  def _bh(self): 
    n = self.nx -1
    dp = ones(n)*2.
    dl = ones(n)*-1.
    ds = dl
  
    dh2 = 1./(self.dx*self.dx)
    D = [dp,dl,ds]
    
    bh = dia_matrix( (D,array([0,-1,1])),shape=(n,n))
    I  = dia_matrix( (-1.*dl,array([0])),shape=(n,n))
    
    A = I + self.a*self.dt*dh2*bh
    return A
    
  def _rhs(self,it):
    b = np.zeros(self.nx-1)
    for ix in range(1,self.nx-1,1):
      b[ix-1] = self.f[it,ix]*self.dt+self.u[it-1,ix]
    return b
  ######## end of implicit Euler #############








class CrankNicolson:
  '''
  '''
  def __init__(self,f,v,g1,g2,st,sx,a=0.1,verb=True):
    ''' 
    init object
    '''
    self.nt = st[2] 
    self.ot = st[0]
    self.dt = st[1]

    self.nx = sx[2]
    self.ox = sx[0]
    self.dx = sx[1]
    self.a  = a
    self.u = np.zeros((self.nt,self.nx))
    self.f = f 
    self.v = v(sx)    
    self.g1 = g1(st)
    self.g2 = g2(st)

    self.verb = verb
    self.alpha = self.a*self.dt/(self.dx*self.dx)


  def solve(self):
    self.u[0,:] = self.v
    A = self._getA()
    for it in range(1,self.nt,1):
      b = self._rhs(it)       
      x = linalg.spsolve(A,b) # linear solver for sparse A

      for ix in range(1,self.nx-1,1):
        self.u[it,ix] = x[ix-1]
      self.u[it,0] = self.g1[it]
      self.u[it,self.nx-1] = self.g2[it]
    return self.u

  ######## Private functions within the object #########
  def _getA(self,lhs=True):
    n = self.nx-2 
    bh = self._getBh() 
    I  = dia_matrix( (ones(n),array([0])),shape=(n,n))
    if lhs:
      A = I + 0.5*self.alpha*bh
    else:
      A = I - 0.5*self.alpha*bh
    return A

  def _getBh(self):
    n = self.nx-2 
    dp = ones(n)*2.
    dl = ones(n)*-1.
    ds = dl
    D = [dp,dl,ds]
    bh = dia_matrix( (D,array([0,-1,1])),shape=(n,n))
    return bh 
  
  def _rhs(self,it):
    um1 = self.u[it-1,1:self.nx-1]
    A = self._getA(False)
    sx = (self.ox,self.dx,self.nx)

    t = (it-0.5)*self.dt+self.ot
    Ftmp = self.dt*self.f(sx,t)

    F = np.zeros(self.nx-1)
    F = Ftmp[1:self.nx-1]

    b = A*um1 +F
     
    b[0] += 0.5*self.alpha*(self.g1[it]+self.g1[it-1])
    b[len(b)-1] += 0.5*self.alpha*( self.g2[it]+self.g2[it-1])
    return b
  ######## end of implicit Euler #############  
