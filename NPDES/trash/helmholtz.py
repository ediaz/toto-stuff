"""
2D Helmholtz solver
"""
import scipy.sparse,scipy.sparse.linalg,scipy.special,scipy.signal
from math import *
from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d

############################################################################

def main():
  #goTestSolution()
  #goCompare()
  goMarmousi()
  #goRateOfConvergence()
  #goTestSolution()
def setGlobals(name=None,h=0.01):
  """Sets the global variables."""
  global nx,dx,fx,ny,dy,fy,freq,c,kxsource,kysource,fname,sigma
  fname = name
  if name=='marmousi':
    nx,dx,fx = 500,0.016,0.000
    ny,dy,fy = 150,dx,fx
    kxsource = 281
    kysource = 3
    sigma = 0.04
    freq = 20.0
    c = 1.0
  else:
    nx = int(1.0/h)-1; dx,fx = h,h
    ny,dy,fy = nx,h,h
    kxsource = 1+nx/2
    kysource = 1+ny/2
    sigma = 0.05
    freq = 6.0
    c = 1.0

############################################################################

def goHelmholtz():
  setGlobals()
  u = solveHelmholtz()
  pixels(scipy.real(u))
  pixels(scipy.imag(u))
  return u

def solveHelmholtz(pml=True):
  a = makeLhs(pml=pml)
  b = makeRhs()
  u = scipy.sparse.linalg.spsolve(a,b)
  return u.reshape(ny,nx)

def makeRhs():
  """Makes the RHS source vector."""
  kx,ky = kxsource,kysource
  b = zeros(nx*ny)
  def gaussian(x,y,sigma=0.05):
    x0 = fx+kx*dx
    y0 = fy+ky*dy
    xm = x-x0
    ym = y-y0
    xms = xm*xm
    yms = ym*ym
    s = 0.5/(sigma*sigma)
    return exp(-xms*s-yms*s)
  for ix in range(nx):
    for iy in range(ny):
      x = fx+ix*dx
      y = fy+iy*dy
      i = ix+iy*nx
      b[i] = gaussian(x,y)
  return b

def makeLhs(rx=None,ry=None,pml=True):
  """Makes the LHS impedance matrix."""
  w = 2.0*pi*freq
  k = w/c; ks = k*k
  odx = ody = 1.0/dx
  odxs = odys = odx*odx
  if not pml:
    rx = ones(nx)
    ry = ones(ny)
  elif rx==None and ry==None:
    rx = makePml(nx)
    ry = makePml(ny)

  def makeB(r):
    n = len(r)
    l,m,u = zeros(n),zeros(n),zeros(n)
    for i in range(n):
      ri = r[i  ]
      rm = r[i-1] if i>0   else r[0]
      rp = r[i+1] if i<n-1 else r[n-1]
      #rm = r[i-1] if i>0   else 0.0
      #rp = r[i+1] if i<n-1 else 0.0
      rs = ri*ri
      dr = rp-rm
      if i>0:
        l[i-1] = (rs-0.25*ri*dr)*odxs
      if i<n-1:
        u[i+1] = (rs+0.25*ri*dr)*odxs
      m[i] = 0.5*ks-2.0*rs*odxs
    return scipy.sparse.spdiags([l,m,u],[-1,0,1],n,n)

  bx = makeB(rx)
  by = makeB(ry)
  eyex = identity(nx)
  eyey = identity(ny)
  a = scipy.sparse.kron(eyey,bx)+scipy.sparse.kron(by,eyex)

  if fname=='marmousi':
    m = getMarmousi()
    for ix in range(nx):
      for iy in range(ny):
        i = ix+iy*nx
        rxi = rx[ix]
        ryi = ry[iy]
        rxs = rxi*rxi
        rys = ryi*ryi
        k = w/m[iy,ix]
        ks = k*k
        a[i,i] = ks-2.0*rxs*odxs-2.0*rys*odys

  return a

def makePml(n,s=None):
  """Computes 1/(1+i*sigma(x)/omega)^2."""
  def getWeights():
    smax = 50.0 # max sigma
    l = int(c/(freq*dx)) # pml width in samples
    ls = l*l
    r = zeros(n)
    for i in range(n):
      x = min(i,n-1-i)
      d = x-l
      if d>0.0:
        pass
      else:
        r[i] = smax*d*d/ls
    return r
  w = 2.0*pi*freq
  r = getWeights()
  for i in range(n):
    ri = r[i]
    if s!=None:
      s[i] = ri
    r[i] = 1.0/(1.0+1.0j*ri/w)
  
  fig = pyplot.figure(n)
  ax = fig.add_subplot(111)
  ax.plot(r.real)
  ax.plot(r.imag)

  return r

############################################################################

def goCompare():
  """Compares solution with Green's function to computed solution."""
  setGlobals()
  u = solveHelmholtz()
  z = solveAnalytic()
  pixels(scipy.real(u))
  pixels(scipy.real(z))
  pixels(scipy.imag(u))
  pixels(scipy.imag(z))

def solveAnalytic():
  """Analytic solution for a Gaussian kernel."""

  s = int(3.0*sigma/dx)
  nxp = nx+2*s
  nyp = ny+2*s
  fxp = fx-s*dx
  fyp = fy-s*dy

  def greens():
    """Green's function."""
    kx = kxsource+s
    ky = kysource+s
    x0 = fxp+kx*dx
    y0 = fyp+ky*dy
    w = 2.0*pi*freq
    k = w/c
    g = zeros((nyp,nxp))
    for ix in range(nxp):
      for iy in range(nyp):
        x = fxp+ix*dx
        y = fyp+iy*dy
        xm = x-x0
        ym = y-y0
        r = sqrt(xm*xm+ym*ym)
        g[iy,ix] = -0.25j*scipy.special.hankel1(0,k*r) if r>0.0 else -0.0
    return g

  def kernel():
    """Gaussian kernel."""
    n = 2*s+1
    f = -s*dx
    k = zeros((n,n))
    def gaussian(x,y):
      xs = x*x
      ys = y*y
      s = 0.5/(sigma*sigma)
      return exp(-xs*s-ys*s)
    for i in range(n):
      for j in range(n):
        x = f+i*dx
        y = f+j*dy
        k[j,i] = gaussian(x,y)
    pixels(scipy.real(k))
    return k

  g = greens()
  k = kernel()
  z = scipy.signal.convolve2d(g,k,mode='same')
  y = z[s:nyp-s,s:nxp-s]
  return y

############################################################################

def goTestSolution(h=0.01):
  """Test with a known solution."""

  setGlobals(h=h)

  # Exact solution
  def u(x,y):
    return sin(pi*x)*sin(pi*y)
  def u_x(x,y):
    return pi*cos(pi*x)*sin(pi*y)
  def u_y(x,y):
    return pi*sin(pi*x)*cos(pi*y)
  def u_xx(x,y):
    return -pi*pi*sin(pi*x)*sin(pi*y)
  def u_yy(x,y):
    return -pi*pi*sin(pi*x)*sin(pi*y)

  # PML function and derivatives
  def rx(x):
    return 1*sin(2.0*pi*x)
  def rx_x(x):
    return 1*2.0*pi*cos(2.0*pi*x)
  def ry(y):
    return 1*sin(2.0*pi*y)
  def ry_y(y):
    return 1*2.0*pi*cos(2.0*pi*y)

  def testSolution():
    """Exact solution."""
    z = zeros(nx*ny)
    for ix in range(nx):
      for iy in range(ny):
        x = fx+ix*dx
        y = fy+iy*dy
        i = ix+iy*nx
        z[i] = u(x,y)
    return z

  def testRhs():
    """Makes the RHS vector."""
    w = 2.0*pi*freq
    k = w/c; ks = k*k
    b = zeros(nx*ny)
    for ix in range(nx):
      for iy in range(ny):
        x = fx+ix*dx
        y = fy+iy*dy
        i = ix+iy*nx
        rxi = rx(x)
        ryi = ry(y)
        rxs = rxi*rxi
        rys = ryi*ryi
        drx = rx_x(x)
        dry = ry_y(y)
        b[i] = ks*u(x,y)+rxi*drx*u_x(x,y)+rxs*u_xx(x,y)+\
                         ryi*dry*u_y(x,y)+rys*u_yy(x,y) 
    return b

  def testLhs():
    """Makes the LHS matrix."""
    w = 2.0*pi*freq
    k = w/c; ks = k*k
    odx = ody = 1.0/dx
    odxs = odys = odx*odx

    lx,mx,ux = zeros(nx),zeros(nx),zeros(nx)
    for ix in range(nx):
      x = fx+ix*dx
      rxi = rx(x)
      rxs = rxi*rxi
      drx = rx_x(x)
      if ix>0:
        lx[ix-1] = rxs*odxs-0.5*rxi*drx*odx
      if ix<nx-1:
        ux[ix+1] = rxs*odxs+0.5*rxi*drx*odx
      mx[ix] = 0.5*ks-2.0*rxs*odxs
    bx = scipy.sparse.spdiags([lx,mx,ux],[-1,0,1],nx,nx)

    ly,my,uy = zeros(ny),zeros(ny),zeros(ny)
    for iy in range(ny):
      y = fy+iy*dy
      ryi = rx(y)
      rys = ryi*ryi
      dry = ry_y(y)
      if iy>0:
        ly[iy-1] = rys*odys-0.5*ryi*dry*ody
      if iy<ny-1:
        uy[iy+1] = rys*odys+0.5*ryi*dry*ody
      my[iy] = 0.5*ks-2.0*rys*odys
    by = scipy.sparse.spdiags([ly,my,uy],[-1,0,1],ny,ny)

    eyex = identity(nx)
    eyey = identity(ny)
    return scipy.sparse.kron(eyey,bx)+scipy.sparse.kron(by,eyex)

  def testPml():
    """Makes the test PML functions."""
    rrx = zeros(nx)
    rry = zeros(ny)
    for ix in range(nx):
      x = fx+ix*dx
      rrx[ix] = rx(x)
    for iy in range(ny):
      y = fy+iy*dy
      rry[iy] = ry(y)
    pyplot.figure()
    pyplot.plot(scipy.real(rrx))
    pyplot.figure()
    pyplot.plot(scipy.real(rry))
    return rrx,rry

  rrx,rry = testPml()
  a = testLhs()
  b = testRhs()
  y = scipy.sparse.linalg.spsolve(a,b)
  z = testSolution()
  pixels(scipy.real(reshape(y)))
  pixels(scipy.real(reshape(z)))
  pixels(scipy.real(reshape(y-z)))
  e = scipy.sum(scipy.absolute(y-z))/(nx*ny)
  #print 'e_h =',e
  return e

############################################################################

def goMarmousi():
  """Solution for Marmousi model."""
  setGlobals('marmousi')
  u = solveHelmholtz()
  pixels(scipy.real(u))
  pixels(scipy.imag(u))
  return u

def getMarmousi():
  v = scipy.fromfile('./marmousi.dat')
  v = v.reshape(ny,nx)
  pixels(v)
  return v

############################################################################

def goRateOfConvergence():
  """Estimates rate of convergence."""
  #hs = [0.02,0.01,0.005]
  hs = [0.02,0.01,0.005,0.0025]
  nh = len(hs)
  e = scipy.zeros(nh)
  for i in range(nh):
    e[i] = goTestSolution(hs[i])
    if i>0:
      r = log10(e[i-1]/e[i])/log10(2.0)
      print 'h=%.4f, e=%.4e, r=%.4f'%(dx,e[i],r)
    else:
      print 'h=%.4f, e=%.4e'%(dx,e[i])

############################################################################
# Utilities

fontsize = 16
def pixels(z):
  if fname=='marmousi':
    fig = pyplot.figure(figsize=(12.0,5.0))
  else:
    fig = pyplot.figure(figsize=(12.0,8.0))
  ax = fig.add_subplot(111)
  cax = ax.imshow(z,interpolation='nearest',origin='upper')
  ax.set_xlabel('x',fontsize=fontsize)
  ax.set_ylabel('y',fontsize=fontsize)
  ax.tick_params(labelsize=fontsize)
  cbar =  fig.colorbar(cax,orientation='vertical',shrink=0.58,aspect=10.0)
  cbartick = pyplot.getp(cbar.ax.axes,'yticklabels')
  pyplot.setp(cbartick,fontsize=fontsize)

def zeros(size):
  return scipy.zeros(size,dtype=complex)

def ones(size):
  return scipy.ones(size,dtype=complex)

def identity(size):
  return scipy.sparse.identity(size,dtype=complex)

def reshape(x):
  return x.reshape(ny,nx)

############################################################################
if __name__ == "__main__":
  main()
  pyplot.show()
