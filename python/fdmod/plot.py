import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import numpy as np

plt.rcParams.update({'font.size': 20,'legend.fontsize': 20})

def plotvel(rsfv,fn=1,interp='bilinear',vmin=-1,vmax=-1,cbar=True,fname=None):
  oz,dz,nz = rsfv.o1, rsfv.d1, rsfv.n1
  ox,dx,nx = rsfv.o2, rsfv.d2, rsfv.n2
  vel = rsfv.rsf 
  ratio = ((nx-1)*dx)/((nz-1)*dz)
  
  fig = plt.figure(fn,figsize=(10,5))
  ax  = fig.add_subplot(111)
  ax.set_xlabel('x(%s)'%rsfv.Aunit)
  ax.set_ylabel('z(%s)'%rsfv.Aunit)

  extent = [ox,ox+(nx-1)*dx,oz+(nz-1)*dz,oz]

  if vmin == -1 or vmax == -1:
    vmin = rsfv.rsf.min()
    vmax = rsfv.rsf.max()


  im = ax.imshow(rsfv.rsf.T,extent=extent,interpolation=interp,vmin=vmin,vmax=vmax)
  if cbar:
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    cb = fig.colorbar(im,cax=cax)
    cb.set_label('%s(%s)'%(rsfv.prop,rsfv.Punit))
  if fname:
    plt.savefig(fname,bbox_inches='tight')



def plotwavefield(rsfw,fn=1,s=0.99,fname=None):
  oz,dz,nz = rsfw.o1, rsfw.d1, rsfw.n1
  ox,dx,nx = rsfw.o2, rsfw.d2, rsfw.n2
  rw = rsfw.rsf.real
  ri = rsfw.rsf.imag
  m = abs(rsfw.rsf).max()

  extent = [ox,ox+(nx-1)*dx,oz+(nz-1)*dz,oz]
  fig = plt.figure(fn,figsize=(12,8))
  ax1 = fig.add_subplot(211)
  #ax1.set_xlabel('x(%s)'%rsfw.Aunit)
  ax1.set_ylabel('z(%s)'%rsfw.Aunit)
  ax1.imshow(rw.T,extent=extent,cmap=cm.gist_gray_r,vmin=-s*m,vmax=s*m) 
  ax2 = fig.add_subplot(212)
  ax2.set_xlabel('x(%s)'%rsfw.Aunit)
  ax2.set_ylabel('z(%s)'%rsfw.Aunit)
  ax2.imshow(ri.T,extent=extent,cmap=cm.gist_gray_r,vmin=-s*m,vmax=s*m)
  if fname:
    plt.savefig(fname,bbox_inches='tight')


def plot1d(rsfv,fn=1,fname=None):
  o,d,n = rsfv.o1,rsfv.d1,rsfv.n1
  x = np.linspace(o,(n-1)*d+o,n)
  fig = plt.figure(fn)
  ax  = fig.add_subplot(111)
  ax.plot(x,rsfv.rsf)  
  ax.set_xlabel('%s(%s)'%(rsfv.Aname,rsfv.Aunit))
  if fname:
    plt.savefig(fname,bbox_inches='tight')


def show():
  plt.show()
