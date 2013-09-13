import numpy as np
from scipy import *




class Waveequation2d:

  def __init__(self,v,M,N,R,c=2.,jsnap=10):
    self.nt = M  # number of grid points
    self.ot = 0.
    self.dt = 1./(M-1)
    self.ntsnap = M/jsnap+1
    self.dtsnap = self.dt*jsnap

    self.ox = 0.
    self.nx = N
    self.dx = 1./(N-1)

    self.ih2 = 1./(self.dx*self.dx)

    self.oy = 0.
    self.ny = self.nx
    self.dy = self.dx

    self.jsnap = jsnap # how many snapshots to avoid for the movie?

    self.c = c
    self.c2 = c*c
    self.v = v
    self.R = R

    self.check_CFL()

    self.um1  = np.zeros((self.nx,self.ny))    
    self.u0   = np.zeros((self.nx,self.ny))    
    self.up1  = np.zeros((self.nx,self.ny)) 
    self.utmp = np.zeros((self.nx,self.ny)) 
    
    self.umovie = np.zeros((self.ntsnap,self.nx,self.ny)) 
    self.get_uo()
    self.get_u1()

  def check_CFL(self):
    if self.c*(self.dt/self.dx)>1:
      print "CFL is violated"

  def solve(self):
    dt2 = self.dt*self.dt

    isnap = 0
    self.umovie[isnap] = self.um1
 
    for it in range(1,self.nt-1,1):
      self.laplacian()
      self.up1 = self.utmp*dt2*self.c2 -self.um1 +2.*self.u0

      # save snapshot to movie
      if ((it+1)%self.jsnap == 0):
        isnap+=1
        self.umovie[isnap] = self.u0
        #print it,self.nt 
      # circulate arrays:
      self.um1 = self.u0
      self.u0 = self.up1 
    return self.umovie



  def get_uo(self):
    x = np.linspace(0.,1.,self.nx)   
    y = np.linspace(0.,1.,self.ny)  

    for ix in range(self.nx):
      for iy in range(self.nx):
        self.um1[ix,iy] = self.v(self.R,x[ix],y[iy])

  
  def get_u1(self):
    c2 = self.c2
    dt = self.dt
    self.u0 = self.um1 +c2*dt*dt*0.5*self.laplacianu0(self.um1)


  
  def laplacian(self):
    for ix in range(1,self.nx-1,1):
      for iy in range(1,self.nx-1,1):
        self.utmp[ix,iy] = self.u0[ix-1,iy] +self.u0[ix+1,iy]+\
                           self.u0[ix,iy-1] +self.u0[ix,iy+1]-\
                           4.*self.u0[ix,iy]
    self.utmp *= self.ih2

  def laplacianu0(self,u):
    uout = np.zeros((self.nx,self.ny))

    for ix in range(1,self.nx-1,1):
      for iy in range(1,self.nx-1,1):
        uout[ix,iy] = u[ix-1,iy] +u[ix+1,iy]+\
                      u[ix,iy-1] +u[ix,iy+1]-\
                           4.*u[ix,iy]
    uout *= self.ih2
    return uout

 
   
