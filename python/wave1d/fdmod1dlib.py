import numpy as np
from scipy import signal

class wave1d:
  def __init__(self,ox,dx,nx,ot,dt,nt,v,rho,nb=0):
    self.ox = ox
    self.dx = dx
    self.nx = nx
    self.nxb = nx+nb*2
    self.nb = nb

    self.ot = ot
    self.dt = dt
    self.nt = nt

    self.v = self.pad(v)
    self.rho = self.pad(rho)
    self.movie = np.zeros((self.nx,nt))

  def pad(self,array):
    out = np.zeros(self.nxb)
    out[0:self.nb] = array[0]
    out[self.nb:self.nx+self.nb]=array
    out[self.nb+self.nx:self.nxb]=array[self.nx-1]
    return out

  def init_source(self,source,sx):
    self.source = source
    self.sx = sx 

  def goforward(self):
    um1  = np.zeros(self.nxb)
    uo   = np.zeros(self.nxb)
    up1  = np.zeros(self.nxb)
    utmp = np.zeros(self.nxb)
    
    rhov2 = self.rho*self.v*self.v*self.dt*self.dt
    buoyancy = np.ones(self.nxb)*np.power(1./self.dx,2)
    
    for ix in range(self.nxb-1): 
      buoyancy[ix] = 2.0/(self.rho[ix]+self.rho[ix+1])
    buoyancy[self.nxb-1] = 1/self.rho[self.nxb-1]
    buoyancy *= np.power(1./self.dx,2)
   
    for it in range(self.nt): 
      self.source_inject(it,uo)
      #
      utmp = self.div_rho_grad(uo, buoyancy)
      #
      up1 = 2.0*uo + rhov2*utmp -um1
      #
      self.movie[:,it] = uo[self.nb:self.nx+self.nb]
      #
      um1 = uo
      uo  = up1

  def div_rho_grad(self,uo,buoyancy):
    u = np.zeros(self.nxb)
    gi = 0.0
    for ix in range(1,self.nxb,1):
      gi       = uo[ix  ] 
      gi      -= uo[ix-1]
      gi      *= buoyancy[ix]
      u[ix-1] -= gi
      u[ix  ]  = gi  
    return -1*u

  def source_inject(self,it,u):
    for isou in range(len(self.sx)):
      ix = int(round((self.sx[isou] - self.ox)/self.dx))
      u[ix+self.nb] += self.source[isou,it]

  def absorb(self,u):
    s = self.nb*2.0
    for ix in range(self.nb):
      f = ix/(np.sqrt(2.)*s)
      u[ix] *= np.exp(-f*f)
      u[self.nx+self.nb+ix] *= np.exp(-f*f)

  def apply(self):
    self.goforward()
    return self.movie 


def ricker(ot,dt,nt,f,t0):
  t = np.linspace(ot,ot+(nt-1)*dt,nt)
  s = 1./(4*f)
  scale = 2./(np.sqrt(3*s)*np.pi**0.25)
  is2 = 1/(2.*s*s)
  t2 = t-t0
  st = -t2*t2*2*is2
  expt = np.exp(-t2*t2*is2)
  return scale*(1+st)*expt

def dgauss(ot,dt,nt,f,t0):
  t = np.linspace(ot,ot+(nt-1)*dt,nt)
  s = 1./(4*f)
  scale = 1/np.sqrt(2*np.pi*s)
  is2 = 1/(2.*s*s)
  t2 = t-t0
  expt = np.exp(-t2*t2*is2)
  return ((-2*(t2)/(2*s*s))*expt)*1/max(((-2*(t2)/(2*s*s))*expt))

