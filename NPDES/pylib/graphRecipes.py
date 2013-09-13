from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np




def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def surfplot(nt,ot,dt,\
             nx,ox,dx,\
             A,fn=1,title=""):
  x = np.arange(ot,ot+nt*dt,dt)
  y = np.arange(ox,ox+nx*dx,dx)

  fig = plt.figure(fn)
  ax  = fig.gca(projection='3d')

  x,y = np.meshgrid(x,y)

  surf = ax.plot_surface(x,y,A, cmap=cm.coolwarm,rstride=1, cstride=1,
                         linewidth=0)

  ax.set_zlim(-1.01, 1.01)
  ax.set_title(r''+title)
  ax.zaxis.set_major_locator(LinearLocator(10))
  ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
  fig.colorbar(surf, shrink=0.5, aspect=5)


def display2d(nt,ot,dt,nx,ox,dx,A,fn=1,title=""):
  x = np.arange(ot,ot+nt*dt,dt)
  y = np.arange(ox,ox+nx*dx,dx)
  mt = x[len(x)-1]  
  mx = y[len(y)-1]  
  x,y = np.meshgrid(x,y)
  

  fig = plt.figure(fn)
  ax = fig.add_subplot(1,1,1)
  ax.imshow(A,cmap=cm.jet,extent=[ot,mt,mx,ox])
  forceAspect(ax,aspect=0.75)
  ax.set_xlabel('x')
  ax.set_ylabel('t')
  ax.set_title(title)

def graph_f(f,nt,ot,dt,fn=1,title=""):
  fig = plt.figure(fn)
  ax = fig.add_subplot(1,1,1)
  x = np.arange(ot,ot+nt*dt,dt)
  plt.plot(x,f)  
  ax.set_title(title)
  ax.set_xlabel('t(s)')


def graph_A(A,nt,ot,dt,j,fn=1,dx=1,title=""):
  fig = plt.figure(fn)
  x =  np.arange(ot,ot+nt*dt,dt)
  n = len(A)
  ax = fig.add_subplot(1,1,1)    
  for it in range(0,n,j): 
    plt.plot(x,A[it],label='t=%3.2gs'%(it*dx))  
    plt.hold(True)
  ax.set_xlabel('x')
  ax.set_title(r''+title)
  plt.legend(loc='upper right')

def show():
  plt.show()


