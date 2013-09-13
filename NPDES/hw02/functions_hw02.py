import numpy as np



def f1(sx,t):
  ox,dx,nx = sx[0],sx[1],sx[2]
  
  x = np.arange(ox,ox+nx*dx,dx)
  u = np.zeros((nx))
  for ix in range(nx):
    u[ix] = np.exp(-2*t)*np.sin(np.pi*x[ix])
  return u

def v1(sx):
  ox,dx,nx = sx[0],sx[1],sx[2]
  x = np.arange(ox,ox+nx*dx,dx)
  u = np.zeros((nx))
  for ix in range(nx):
    u[ix] = x[ix]*(1-x[ix]) 
  return u


def f2(sx,t):
  ox,dx,nx = sx[0],sx[1],sx[2]
  x = np.arange(ox,ox+nx*dx,dx)
  u = np.zeros((nx))
  return u


def v2(sx):
  ox,dx,nx = sx[0],sx[1],sx[2]
  x = np.arange(ox,ox+nx*dx,dx)
  u = np.zeros((nx))
  for ix in range(nx):
    u[ix] = np.sin(2.*np.pi*x[ix]) 
  return u



def uexact(sx,st):
  ox,dx,nx = sx[0],sx[1],sx[2]
  ot,dt,nt = st[0],st[1],st[2]
  
  u = np.zeros((nt,nx))
  a = 1./12.
  for it in range(nt):
    t = ot+dt*it
    xx = -a*(2.*np.pi)*(2.*np.pi)*t
    for ix in range(nx):  
      x = ox+dx*ix
      u[it,ix] = np.sin(2.*np.pi*x)*np.exp(xx)
  return u


def l2err(err):
  err2 = err*err
  N = len(err)
  return np.sqrt(np.sum(err2))/N
