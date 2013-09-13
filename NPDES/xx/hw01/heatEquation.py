import numpy as np
from scipy import *
from scipy.sparse import dia_matrix,linalg

'''
This module provides the Crank-Nicolson method for
solving the 1d heat-diffusion equation.

Crank-Nicolson class:
  solves the 1d heat equation by solving the linear system of 
  equations:
 
  [I + 0.5alpha] Um = [I-0.5alpha]Um-1 dtF^{m-0.5} +0.5alpha C

'''


class CrankNicolson:
  '''
  '''
  def __init__(self,f,v,g1,g2,st,sx,a=0.1,verb=True):
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
    self.g1 = g1(st)
    self.g2 = g2(st)

    self.verb = verb
    self.alpha = self.a*self.dt/(self.dx*self.dx)
  
  # The problem is solved in the following routine
  def solve(self):
    self.u[0,:] = self.v
    A = self._getA()

    for it in range(1,self.nt,1):
      b = self._rhs(it)       
      x = linalg.spsolve(A,b) # linear solver for sparse A

      self.u[it,1:self.nx-1] = x
      self.u[it,0] = self.g1[it]
      self.u[it,self.nx-1] = self.g2[it]

    return self.u




  ######## Private functions within the object #########
  def _getA(self,lhs=True):
    '''
    gets the left hand side of the linear system of 
    equations.
    
    If lhs=False it returns the matrix needed     
    in the right hand side
    '''
    n = self.nx-2 
    bh = self._getBh() 
    I  = dia_matrix( (ones(n),array([0])),shape=(n,n))
    if lhs:
      A = I + 0.5*self.alpha*bh
    else:
      A = I - 0.5*self.alpha*bh
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
