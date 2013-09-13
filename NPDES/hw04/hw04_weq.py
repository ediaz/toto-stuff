import numpy as np
import waveequation as weq
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import os

def surf(u1,obj,title="",figname=None):
  nx = obj.nx
  ny = obj.ny
  X = np.linspace(0.,1.,nx)
  Y = np.linspace(0.,1.,ny)
  X, Y = np.meshgrid(X, Y)
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  surf = ax.plot_surface(X, Y, u1, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
  ax.set_zlim(-1.01, 1.01)
  ax.set_title(title)
  if not figname == None:
    fig.savefig(figname)


def plot_side(u1,u2,fn):
  scale = 1. 
  fig = plt.figure(fn)
  ax1  = fig.add_subplot(1,2,1)
  ax1.imshow(u1,vmin=-scale,vmax=scale)
  ax2  = fig.add_subplot(1,2,2)
  im2= ax2.imshow(u2,vmin=-scale,vmax=scale)
  cbar = fig.colorbar(im2)


def error(u,u1):

  e = np.sqrt((abs(u-u1)*abs(u-u1)).sum())/len(u)
  return e


def v(R,x,y):
  cx = 0.5
  cy = 0.5
  
  d = np.sqrt((x-cx)**2+(y-cy)**2)
  if R >= d:
    out = np.cos(0.5*np.pi*d/R)
    out *= out
  else:
    out = 0.
  return out





#######################################################################
# Here I create the movie frames for diffusion example
N = 1601
M = (N-1)/4+1
c = 2.0
jsnap = 5


 
outfile1 = 'solution-N%03d-C%3.1f.npy'%(N,c)
#outfile2 = 'exact-N%03d-C%3.1f.npy'%(N,c)
modeling = 'weq/'
folder  = 'solutions/'

if not os.path.exists(folder+modeling+outfile1):
  weq2d= weq.Waveequation2d(v,N,M,0.3,c=2.,jsnap=jsnap)
  u = weq2d.solve()
  np.save(folder+modeling+outfile1,u)
  movieframesdir = folder+modeling+'movieN%03d-C%3.1f'%(N,c)
  os.makedirs(movieframesdir)
  for it in range(0,len(u),10):
    figname = movieframesdir+'/frame-it%02d.png'%it
    tit = "t=%5.3f"%(it*weq2d.dt*weq2d.jsnap)
    surf(u[it],weq2d,title=tit,figname=figname)  


#
nk = 4
#p = ['']*nk
#L2 = 1.

ufine = np.load(folder+modeling+outfile1)
Nf = len(ufine)
Mf = len(ufine[0])
#
#
uex = ufine[:,Mf/2,Mf/2]

e = np.zeros(nk)
e[1] = 1.
k=0
for n in [25,51,101,201]:
  N = 1601 
  M = n
  Ns = N/jsnap+1
  Nr = (Nf-1)/(Ns-1)
  Mr = (Mf-1)/(M-1)


  weq2d= weq.Waveequation2d(v,N,M,0.3,c=2.,jsnap=jsnap)
  umov = weq2d.solve()
  u = umov[:,M/2,M/2] 
  print u.shape
 
  e[k] = error(u,uex) 
  print e[k]
  k+=1

k=0
for n in [25,51,101,201]:
  N = n 
  M = (N-1)/4+1
  dt = 1./(N-1)
  print n,1./(n-1),e[k],np.log(e[k-1]/e[k])/np.log(2)
  k+=1



#  N = int(.1/h)+1
#  CN2d = heq.CrankNicolson2d(g,N,c=c,jsnap=5)
#  u = CN2d.solve()
#  uex = CN2d.uexact() 
#  e = error(u,uex)
#  r = np.log(L2/e)/np.log(2)
#  L2 = e
#  print '%d& %4d& %10.8f& %10.8g& %10.8g'%(k,N,h,e,r)
#
#for k in range(nk):
#  print p[k]
