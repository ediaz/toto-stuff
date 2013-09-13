import numpy as np



def uexact(sx,st):
  ox,dx,nx = sx[0],sx[1],sx[2]
  ot,dt,nt = st[0],st[1],st[2]
  
  x = np.arange(ox,ox+nx*dx,dx)
  t = np.arange(ot,ot+nt*dt,dt)
  u = np.zeros((nt,nx))
  a = 0.1
  for it in range(nt):
    for ix in range(nx):
      u[it,ix] = (np.sin(4.*x[ix])+np.cos(2.*x[ix]))*\
                 (np.cos(t[it])+np.sin(t[it]))

  return u

def f(sx,t):
  ox,dx,nx = sx[0],sx[1],sx[2]
  
  x = np.arange(ox,ox+nx*dx,dx)
  u = np.zeros((nx))
  a=0.1
  for ix in range(nx):
    cost = np.cos(t)
    sint = np.sin(t)
    sin4x = np.sin(4.*x[ix])
    cos2x = np.cos(2.*x[ix])
    u[ix] = cost*((1.+16.*a)*sin4x +(1.+4.*a)*cos2x)+\
               sint*((4.*a-1)*cos2x+(16.*a-1.)*sin4x)
  return u

def v(sx):
  ox,dx,nx = sx[0],sx[1],sx[2]
  x = np.arange(ox,ox+nx*dx,dx)
  u = np.zeros((nx))
  for ix in range(nx):
    u[ix] = (np.sin(4.*x[ix])+np.cos(2.*x[ix]))
  return u
  
def g1(st):
  ot,dt,nt = st[0],st[1],st[2]

  t = np.arange(ot,ot+nt*dt,dt)
  u = np.zeros((nt))
  x = 2. 
  for it in range(nt):
    u[it] = (np.sin(4.*x)+np.cos(2.*x))*\
                 (np.cos(t[it])+np.sin(t[it]))
  return u
    
def g2(st):
  ot,dt,nt = st[0],st[1],st[2]
  t = np.arange(ot,ot+nt*dt,dt)
  u = np.zeros((nt))
  x = 3. 
  for it in range(nt):
    u[it] = (np.sin(4*x)+np.cos(2*x))*\
                 (np.cos(t[it])+np.sin(t[it]))
  return u
    
 
def compute_error(u,sx,st):
  a = uexact(sx,st)
  b = abs(a-u)  
  return abs(b[sx[2]]).max()


