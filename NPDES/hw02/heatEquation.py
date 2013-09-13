import numpy as np
from scipy import *
from scipy.sparse import dia_matrix,linalg


class t3level:
  '''
  here I solve:

  3U^m_n -4U^m-1_n +U^m-1_n     U^m_n+1 -2U^m_n +U^m_n-1
  -------------------------- -a ------------------------ = f^m_n
            2dt                          h^2

  for m>= 2

  For m = 1 we have the Backward Euler solution

  '''
  def __init__(self,f,v,st,sx,a=0.1,verb=True):
    ''' 
    init object:
     
    user pass pointers to f,v,g1,g2
    sx,st contain the discretization geometry
    
    nt,nx is the number of gridpoints rather than the
          number of intervals as is the notation of 
          the class
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

    self.verb = verb
    self.alpha = self.a*self.dt/(self.dx*self.dx)
  
  # The problem is solved in the following routine
  def solve(self):
    self.u[0,:] = self.v
    A1 = self._getA(init=True)
    A  = self._getA(init=False)

    b = self._rhs(1)     

    self.u[1,1:self.nx-1] = linalg.spsolve(A1,b)

    for it in range(2,self.nt,1):
      b = self._rhs(it)       
      self.u[it,1:self.nx-1] = linalg.spsolve(A,b) 
    return self.u




  ######## Private functions within the object #########
  def _getA(self,init=True):
    '''
    gets the left hand side of the linear system of 
    equations for m=1.
    '''
    n = self.nx-2 
    bh = self._getBh() 
    I  = dia_matrix( (ones(n),array([0])),shape=(n,n))
    if init:
      A = I + self.alpha*bh
    else:
      A = 3*I + 2.*self.alpha*bh
    return A

  def _getBh(self):
    '''
    Simply returns the second derivative operator 
    '''
    n = self.nx-2 
    dp = ones(n)*2.
    dl = ones(n)*-1.
    ds = dl
    D = [dp,dl,ds]
    bh = dia_matrix( (D,array([0,-1,1])),shape=(n,n))
    return bh 
  
  def _rhs(self,it):
    '''

    '''
    sx = (self.ox,self.dx,self.nx)
    t = it*self.dt +self.ot

    if it<2:
      b = self.dt*self.f(sx,t)+self.u[it-1,:]
    else:
      b = 2*self.dt*self.f(sx,t)-self.u[it-2,:]+4*self.u[it-1,:]
    return b[1:self.nx-1]
  ######## end of implicit Euler ############# 




class implicit:
  '''
    U^m_n - U^m-1_n              U^m_n+1 -2U^m_n +U^m_n-1
  -------------------------- -a ------------------------ = f^m_n
            2dt                          h^2

  '''
  def __init__(self,f,v,st,sx,a=0.1,verb=True):
    ''' 
    init object:
     
    user pass pointers to f,v,g1,g2
    sx,st contain the discretization geometry
    
    nt,nx is the number of gridpoints rather than the
          number of intervals as is the notation of 
          the class
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

    self.verb = verb
    self.alpha = self.a*self.dt/(self.dx*self.dx)
  
  # The problem is solved in the following routine
  def solve(self):
    self.u[0,:] = self.v
    A = self._getA(init=True)
    for it in range(1,self.nt,1):
      b = self._rhs(it)       
      self.u[it,1:self.nx-1] = linalg.spsolve(A,b) 
    return self.u




  ######## Private functions within the object #########
  def _getA(self,init=True):
    '''
    gets the left hand side of the linear system of 
    equations for m=1.
    '''
    n = self.nx-2 
    bh = self._getBh() 
    I  = dia_matrix( (ones(n),array([0])),shape=(n,n))
    if init:
      A = I + self.alpha*bh
    else:
      A = 3*I + 2.*self.alpha*bh
    return A

  def _getBh(self):
    '''
    Simply returns the second derivative operator 
    '''
    n = self.nx-2 
    dp = ones(n)*2.
    dl = ones(n)*-1.
    ds = dl
    D = [dp,dl,ds]
    bh = dia_matrix( (D,array([0,-1,1])),shape=(n,n))
    return bh 
  
  def _rhs(self,it):
    '''

    '''
    sx = (self.ox,self.dx,self.nx)
    t = it*self.dt +self.ot

    b = self.dt*self.f(sx,t)+self.u[it-1,:]
    return b[1:self.nx-1]
  ######## end of implicit Euler #############   
