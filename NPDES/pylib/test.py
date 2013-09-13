import numpy as np
import heatEquation as heq
from graphRecipes import *

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
    for ix in range(nx):
      x = ox +dx*ix
      t = ot +dt*it
      a[it,ix] = np.sin(np.pi*x)*np.exp(-2.*t)
  return a 

def v1(nx):
  a = np.zeros(nx)
  L = 1.0
  ox,dx,nx = ganesh2Toto(L,nx)
  for ix in range(nx):
    x = ox +dx*ix
    if (x<L/2):    
      a[ix] = x
    else:
      a[ix] = 1.0 -x
  return a

    

M = 501
N = 11
ot,dt,nt = ganesh2Toto(1.,M)
ox,dx,nx = ganesh2Toto(1.,N)


# Exercise sheet 2
a1 = heq.explicitEuler(f1,v1,\
                       nt,ot,dt,\
                       nx,ox,dx,\
                       a=1.0)
sol1 = a1.solve()


display2d(nt,ot,dt,nx,ox,dx,sol1,1)
display2d(nt,ot,dt,nx,ox,dx,f1(nx,nt),2)

graph_f(v1(nx),nx,ox,dx,3)




a2 = heq.implicitEuler(f1,v1,\
                       nt,ot,dt,\
                       nx,ox,dx,\
                       a=1.0)

sol2 = a2.solve()
display2d(nt,ot,dt,nx,ox,dx,sol2,4)

show()
