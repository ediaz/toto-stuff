import numpy as np
import heat2dequation as heq
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import os

def surf(u1,obj,figname):
  nx = obj.nx
  ny = obj.ny
  X = np.linspace(0.,1.,nx)
  Y = np.linspace(0.,1.,ny)
  X, Y = np.meshgrid(X, Y)
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  surf = ax.plot_surface(X, Y, u1, rstride=10, cstride=10, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
  ax.set_zlim(-1.01, 1.01)
  fig.savefig(figname)


def plot_side(u1,u2,fn):
  scale = 1. 
  fig = plt.figure(fn)
  ax1  = fig.add_subplot(1,2,1)
  ax1.imshow(u1,vmin=-scale,vmax=scale)
  ax2  = fig.add_subplot(1,2,2)
  ax2.imshow(u2,vmin=-scale,vmax=scale)

def error(u,u1):
  n3 = len(u)
  n2 = len(u[0])
  n1 = len(u[0,0])

  e = (abs(u-u1)).max()
  return e

def g(c,x,y,t):
  return np.sin(np.pi*c*x)*np.sin(np.pi*c*y)*np.exp(-c*c*np.pi*np.pi*t)


#######################################################################
# Here I create the movie frames for diffusion example
N = 101
c = 2.5
outfile1 = 'solution-N%03d-C%3.1f.npy'%(N,c)
outfile2 = 'exact-N%03d-C%3.1f.npy'%(N,c)
modeling = 'heq/'
folder  = 'solutions/'

if not os.path.exists(folder+modeling+outfile1):
  CN2d = heq.CrankNicolson2d(g,N,c=c,jsnap=5)
  u = CN2d.solve()
  uex = CN2d.uexact() 
  np.save(folder+modeling+outfile1,u)
  np.save(folder+modeling+outfile2,uex)
  movieframesdir = folder+modeling+'movieN%03d-C%3.1f'%(N,c)
  os.makedirs(movieframesdir)
  for it in range(0,len(u),10):
    figname = movieframesdir+'/frame-it%02d.png'%it
    figname2 = movieframesdir+'/frameex-it%02d.png'%it
    surf(u[it],CN2d,figname)  
    surf(uex[it],CN2d,figname2)  


nk = 5
p = ['']*nk
L2 = 1.
for k in range(1,nk,1):
  h = 0.02/2**k
  N = int(.1/h)+1
  CN2d = heq.CrankNicolson2d(g,N,c=c,jsnap=5)
  u = CN2d.solve()
  uex = CN2d.uexact() 
  e = error(u,uex)
  r = np.log(L2/e)/np.log(2)
  L2 = e
  print '%d& %4d& %10.8f& %10.8g& %10.8g'%(k,N,h,e,r)

