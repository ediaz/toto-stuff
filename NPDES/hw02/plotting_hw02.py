import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np



plt.rcParams.update({'font.size': 20,'legend.fontsize': 20})

def forceAspect(ax,aspect=1):
  im = ax.get_images()
  extent =  im[0].get_extent()
  ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)




def display3d_3(nt,ot,dt,nx,ox,dx,A1,A2,A3,fn=1,title="",save=None):
  from mpl_toolkits.mplot3d import Axes3D
  from matplotlib import cm
  from matplotlib.ticker import LinearLocator, FormatStrFormatter
  import matplotlib.pyplot as plt

  mt = ot+(nt)*dt  
  mx = ox+(nx)*dx
  
  x=np.zeros(nt)
  y=np.zeros(nx)

  for it in range(nt):
    x[it] = ot+it*dt

  for ix in range(nx):
    y[ix] = ox+ix*dx

  x,y = np.meshgrid(x,y)
    
  fig = plt.figure(fn,figsize=(15,5))
  ax = fig.gca(projection='3d')

  print x.shape,y.shape,A1.shape

  surf = ax.plot_surface(x, y, A1, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)


def display2d(nt,ot,dt,nx,ox,dx,A1,fn=1,title="",save=None):
  mt = ot+(nt-1)*dt  
  mx = ox+(nx-1)*dx
  x = np.arange(ot,mt,dt)
  y = np.arange(ox,mx,dx)
  x,y = np.meshgrid(x,y)

  maxval=(np.array([A1])).max()
  minval=(np.array([A1])).min()

 
  fig = plt.figure(fn,figsize=(5,5))
  ax1 = fig.add_subplot(1,1,1)
  cax1 = ax1.imshow(A1,cmap=cm.jet,extent=[ox,mx,mt,ot])

  # left plot: solution
  forceAspect(ax1,aspect=0.75)
  ax1.set_xlabel('x')
  ax1.set_ylabel('t')
  ax1.set_xticks([ox,0.5*(ox+mx),mx])



def display2d_3(nt,ot,dt,nx,ox,dx,A1,A2,A3,fn=1,title="",save=None):
  mt = ot+(nt-1)*dt  
  mx = ox+(nx-1)*dx
  x = np.arange(ot,mt,dt)
  y = np.arange(ox,mx,dx)
  x,y = np.meshgrid(x,y)

  maxval=(np.array([A1,A2])).max()
  minval=(np.array([A1,A2])).min()

 
  fig = plt.figure(fn,figsize=(15,5))
  ax1 = fig.add_subplot(1,3,1)
  cax1 = ax1.imshow(A1,cmap=cm.jet,extent=[ox,mx,mt,ot],vmin=-maxval,vmax=maxval)

  # left plot: solution
  forceAspect(ax1,aspect=0.75)
  ax1.set_xlabel('x')
  ax1.set_ylabel('t')
  ax1.set_xticks([ox,0.5*(ox+mx),mx])

  cbar = fig.colorbar(cax1,ticks=[maxval,0.5*(maxval+minval),minval],\
                      format='%1.2g',shrink=0.92)

  # center plot: exact solution
  ax2 = fig.add_subplot(1,3,2)
  cax2 = ax2.imshow(A2,cmap=cm.jet,extent=[ox,mx,mt,ot],vmin=-maxval,vmax=maxval)

  forceAspect(ax2,aspect=0.75)
  ax2.set_yticks([])
  ax2.set_xlabel('x')
  cbar = fig.colorbar(cax2,ticks=[maxval,0.5*(maxval+minval),minval],\
                      format='%1.2g',shrink=0.92)
  ax2.set_xticks([ox,0.5*(ox+mx),mx])

  # third plot: difference
  maxd = A3.max()
  mind = A3.min()
  ax3 = fig.add_subplot(1,3,3)
  cax3 = ax3.imshow(A3,cmap=cm.jet,extent=[ox,mx,mt,ot])
  forceAspect(ax3,aspect=0.75)
  ax3.set_yticks([])
  ax3.set_xlabel('x')
  cbar = fig.colorbar(cax3,ticks=[mind,0.5*(maxd+mind),maxd],\
                      format='%3.1g',shrink=0.92)
  ax3.set_xticks([ox,0.5*(ox+mx),mx])

  if not save==None:
    fig.savefig(save)



def display2d_2(nt,ot,dt,nx,ox,dx,A1,A2,fn=1,title="",save=None):
  mt = ot+(nt-1)*dt  
  mx = ox+(nx-1)*dx
  x = np.arange(ot,mt,dt)
  y = np.arange(ox,mx,dx)
  x,y = np.meshgrid(x,y)

  maxval=(np.array([A1,A2])).max()
  minval=(np.array([A1,A2])).min()

 
  fig = plt.figure(fn,figsize=(10,5))
  ax1 = fig.add_subplot(1,2,1)
  cax1 = ax1.imshow(A1,cmap=cm.jet,extent=[ox,mx,mt,ot],vmin=-maxval,vmax=maxval)

  # left plot: solution
  forceAspect(ax1,aspect=0.75)
  ax1.set_xlabel('x')
  ax1.set_ylabel('t')
  ax1.set_xticks([ox,0.5*(ox+mx),mx])

  cbar = fig.colorbar(cax1,ticks=[maxval,0.5*(maxval+minval),minval],\
                      format='%1.2g',shrink=0.92)

  # center plot: exact solution
  ax2 = fig.add_subplot(1,2,2)
  cax2 = ax2.imshow(A2,cmap=cm.jet,extent=[ox,mx,mt,ot],vmin=-maxval,vmax=maxval)

  forceAspect(ax2,aspect=0.75)
  ax2.set_yticks([])
  ax2.set_xlabel('x')
  cbar = fig.colorbar(cax2,ticks=[maxval,0.5*(maxval+minval),minval],\
                      format='%1.2g',shrink=0.92)
  ax2.set_xticks([ox,0.5*(ox+mx),mx])

  if not save==None:
    fig.savefig(save)






def error_plot(xh,eh,fn=1,save=None):
  fig = plt.figure(fn,figsize=(8,8))
  ax1 = plt.subplot(111)
  ax1.set_xscale("log")
  ax1.set_yscale("log")
  ax1.set_xlim(5e-3, 2e-1)
  ax1.set_ylim(4e-5, 5e-2)
  ax1.set_aspect(0.5)
  ax1.set_title("")
  ax1.set_xlabel("h")
  ax1.set_ylabel("e")
  ax1.plot(xh,eh, "o-")
  plt.draw()
  if not save==None:
    fig.savefig(save)

def show():
  plt.show()

