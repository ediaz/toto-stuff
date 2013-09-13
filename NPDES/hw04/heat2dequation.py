import numpy as np
from scipy import *
from scipy.sparse import dia_matrix,linalg,kron







class CrankNicolson2d:

  def __init__(self,g,N,c=2.5,a=0.5,jsnap=10):
    self.nt = N  # number of grid points
    self.ot = 0.
    self.dt = .1/(N-1)
    self.ntsnap = N/jsnap+1
    self.dtsnap = self.dt*jsnap

    self.ox = 0.
    self.nx = int(1./self.dt)+1
    self.dx = self.dt

    self.oy = 0.
    self.ny = self.nx
    self.dy = self.dx

    self.jsnap = jsnap # how many snapshots to avoid for the movie?
    self.n = self.nx-2 # number of unknowns within the domain
    self.g = g
    self.a = a # conductivity
    self.c = c
    self.alpha =  self.a*self.dt/(self.dx*self.dx)

    self.um1  = np.zeros(self.n*self.n)    
    self.u0   = np.zeros(self.n*self.n)    
    self.utmp = np.zeros(self.n*self.n)    
    self.umovie = np.zeros((self.ntsnap,self.nx,self.nx)) 
    self.umovex = np.zeros((self.ntsnap,self.nx,self.nx)) 
    self._fillbounds()

    print 'nt = %d   nx = %d  ny = %d'%(self.nt,self.nx,self.ny)


  def solve(self):
    self._fillt0()
    A = self._getA()
    self._Arhs()
    u = np.reshape(self.um1,(self.n,self.n))
    isnap = 0
    self.umovie[isnap,1:self.nx-1,1:self.ny-1] = u
    for it in range(1,self.nt,1):

      b = self._rhs(it)       
      self.u0 = linalg.spsolve(A,b) # linear solver for sparse A
      if ((it)%self.jsnap == 0):
        print '%d/%d'%(it,self.nt)
        isnap+=1
        u = np.reshape(self.u0,(self.n,self.n))    
        self.umovie[isnap,1:self.nx-1,1:self.ny-1] = u
        
      self.um1 = self.u0  # swap pointers 
      self.u0 = np.zeros(self.n*self.n)
    return self.umovie



  # private methods:

  def _fillt0(self):
    d = self.dx
    n = self.n
    y = x = np.linspace(d,1.-d,n)
    self.um1 = np.reshape(self.um1,(n,n))
    for ix in range(n):  
      for iy in range(n):  
        self.um1[ix,iy] = self.g(self.c,x[ix],y[iy],0.) 
    self.um1 = np.reshape(self.um1,(n*n))
   

  def _getA(self,lhs=True):
    '''
    gets the left hand side of the linear system of 
    equations.
    
    If lhs=False it returns the matrix needed     
    in the right hand side
    '''
    n = self.n
    n2 = n*n
    bh = self._laplace() 
    I  = dia_matrix((ones(n2),array([0])),shape=(n2,n2))
    if lhs:
      A = I + 0.5*self.alpha*bh
    else:
      A = I - 0.5*self.alpha*bh
    return A



  def _laplace(self):
    '''
    construct Laplace operator as a matrix
    based on user input filter.
    Actually this matrix can be thought
    as a convolution operator:
    f(x,z)*U(x,z)
    '''
    f=[-1.,2.,-1.]
    nx = nz = self.n
    nf = len(f)
    nonzero = np.ones((nf,nx))
    for i in range(nf):
      nonzero[i] *=f[i]
    offsets = array(range(nf/2,-nf/2 ,-1))
  
    m1 = dia_matrix((nonzero,offsets),shape=(nx,nx))
    m2 = identity(nz)
    k1 = kron(m1,m2)
    nonzero = np.ones((nf,nz))
    for i in range(nf):
      nonzero[i,:] *=f[i]
    m1 = dia_matrix((nonzero,offsets),shape=(nz,nz))
    m2 = identity(nx)
    k2 = kron(m2,m1)   
    return k1+ k2


  def _gett(self,it):
    return it*self.dt +self.ot

  def _getx(self,ix):
    return ix*self.dx +self.ox

  def _Arhs(self):
    self.Arhs = self._getA(False)

  def _rhs(self,it):
    '''

    '''
    um1 = self.um1 
    b = self.Arhs*um1 
    t0 = self._gett(it) 
    tm1 = self._gett(it-1) 
    b += 0.5*self.alpha*(self._g1(t0)+self._g1(tm1))
    return b

  def _fillbounds(self):
    '''
    this method introduces the boundary data
    '''
    x = np.linspace(0.,1.,self.nx)
    y = np.linspace(0.,1.,self.nx)

    for it in range(len(self.umovie)):
      t = self.ot +it*self.dt*self.jsnap 
      for iy in range(self.ny):
        yy = y[iy]
        self.umovie[it,0,iy] = self.g(self.c,0,yy,t)
        self.umovie[it,self.nx-1,iy] = self.g(self.c,x[self.nx-1],yy,t)

      for ix in range(self.nx):
        xx = x[ix]
        self.umovie[it,ix,0] = self.g(self.c,xx,0,t)
        self.umovie[it,ix,self.ny-1] = self.g(self.c,xx,y[self.ny-1],t)

  
  def _g1(self,t):
    n = self.n
    c = self.c
    u = np.zeros((n,n))
    ux = np.zeros(n)
    uy = np.zeros(n)
    d = 1./(n+1)
    y = x = np.linspace(d,1.-d,n)
    ux0 = self.g(self.c,x,0.,t) 
    ux1 = self.g(self.c,x,1.,t) 
    u0y = self.g(self.c,0.,y,t) 
    u1y = self.g(self.c,1.,y,t)
 
    for ix in range(n):
      u[ix,0  ] = ux0[ix]
      u[ix,n-1] = ux1[ix]

    for iy in range(n):
      u[n-1,iy] += u1y[iy]
      u[0  ,iy] += u0y[iy]
    return np.reshape(u,n**2)


  def uexact(self):
    jsnap = self.jsnap
    x = np.linspace(0.,1.,self.nx)   
    y = np.linspace(0.,1.,self.ny)   
    t = np.linspace(0.,.1,self.ntsnap)   
 
    x,y = np.meshgrid(x,y) 
    umov = np.zeros((self.ntsnap,self.nx,self.ny))
    for it in range(self.ntsnap):
      umov[it] = self.g(self.c,x,y,t[it])
    return umov
  
   
