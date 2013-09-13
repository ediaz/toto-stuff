import numpy as np

def ganesh2Toto(L,N):
  d = float(L/(N-1))
  o = 0.0
  return o,d,N

def f1(nt,nx):
  L = 1.0 
  ox,dx,nx = ganesh2Toto(L,nx)
  ot,dt,nt = ganesh2Toto(L,nt)
  a = np.zeros([nt,nx])
  for it in range(nt):
    t = ot +dt*it
    for ix in range(nx):
      x = ox +dx*ix
      a[it,ix] = np.sin(np.pi*x)*np.exp(-2.*t)
  return a 


def f2(nt,nx):
  L = 1.0 
  ox,dx,nx = ganesh2Toto(L,nx)
  ot,dt,nt = ganesh2Toto(L,nt)
  a = np.zeros([nt,nx])
  for it in range(nt):
    t = ot +dt*it
    for ix in range(nx):
      x = ox +dx*ix
      a[it,ix] = np.sin(np.pi*x)*np.exp(t) 
  return a 

def f3(nt,nx):
  L = 1.0 
  ox,dx,nx = ganesh2Toto(L,nx)
  ot,dt,nt = ganesh2Toto(L,nt)
  a = np.zeros([nt,nx])
  for it in range(nt):
    t = ot +dt*it
    for ix in range(nx):
      x = ox +dx*ix
      a[it,ix] = 64*(x*(1.-x)*t)-4*x
  return a 



def v1(nx):
  a = np.zeros(nx)
  L = 1.0
  ox,dx,nx = ganesh2Toto(L,nx)
  for ix in range(nx):
    x = ox +dx*ix
    if (x<L/2):    
      a[ix] = 2*x
    else:
      a[ix] = 2. -2*x
  return a

def v2(nx):
  a = np.zeros(nx)
  L = 1.0
  ox,dx,nx = ganesh2Toto(L,nx)
  for ix in range(nx):
    x = ox +dx*ix
    a[ix] = np.sin(2*np.pi*x) 
  return a

def v3(nx):
  a = np.zeros(nx)
  L = 1.0
  ox,dx,nx = ganesh2Toto(L,nx)
  for ix in range(nx):
    x = ox +dx*ix
    a[ix] = x*x -x
  return a

###### Excercise 3 ################

def u_exact(nt,nx):
  L = 1.0 
  ox,dx,nx = ganesh2Toto(L,nx)
  ot,dt,nt = ganesh2Toto(L,nt)
  a = np.zeros([nt,nx])
  for it in range(nt):
    t = ot +dt*it
    for ix in range(nx):
      x = ox +dx*ix
      a[it,ix] = (np.sin((np.pi*x)/2.)-x)*np.exp(-(np.pi*np.pi*t/4.))+\
                 0.5*np.sin(2*np.pi*x)*np.exp(-4*np.pi*np.pi*t)
  return a 


def f_exact(nt,nx):
  L = 1.0 
  ox,dx,nx = ganesh2Toto(L,nx)
  ot,dt,nt = ganesh2Toto(L,nt)
  a = np.zeros([nt,nx])
  for it in range(nt):
    t = ot +dt*it
    for ix in range(nx):
      x = ox +dx*ix
      a[it,ix] = x*np.pi*np.pi*0.25*np.exp(-np.pi*np.pi*t/4.)
  return a   


def v_exact(nx):
  a = np.zeros(nx)
  L = 1.0
  ox,dx,nx = ganesh2Toto(L,nx)
  t=0.0 
  for ix in range(nx):
    x = ox +dx*ix
    a[ix] = np.sin((np.pi*x)/2.)-x+0.5*np.sin(2.*np.pi*x)
  return a


def inf_norm(A,ax=None):
  return abs(A).max(axis=ax)


