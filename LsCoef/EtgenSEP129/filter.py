import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

plt.rcParams.update({'font.size': 20,'legend.fontsize': 20})


def getFourierAxis(d,n):
  delta = np.pi/(d*(n-1))
  maxF  = delta*(n-1)
  axis = np.arange(0,maxF+delta,delta)
  return axis


class LSfilt:
  def __init__(self,dt,dx,dz,v,nx,nz):
    self.N = 101
    self.dt = dt
    self.dx = dx
    self.dz = dz
    self.nx = nx
    self.nz = nz
    self.v  = v

  def setW(self,eps):
    '''
    set John's suggested weigthing function:  
  
     w(|k|) =  1/|k|^/{1+\epsilon}
     0<eps<=0.5
    '''
    self.w = 1.0/np.power(self.absk(),1.0+eps)

  def of(self,filtx,filtz):
    vphase = self.vphase(filtx,filtz)
    res = np.power(vphase-np.ones((self.N,self.N))*self.v,2)
    return (self.w*res).sum()


  def hessian(self,inp):
    pass

  def vphase(self,filtx,filtz):
    '''
    implement equation 5 from report
    '''
    v = self.v 
    dt = self.dt
    dx = self.dx
    dz = self.dz 
    kk = self.absk()
    vph = np.arccos((v*v*dt*dt/(2.0*dx*dx))*self.Fx(filtx)+\
                    (v*v*dt*dt/(2.0*dz*dz))*self.Fz(filtz)+\
                    np.ones((self.N,self.N)))/(dt*kk)
    return vph

  def Fx(self,filtx):
    nkx = self.N
    dkx = np.pi/(self.dx*(self.N-1))
    xx = np.zeros((self.N,self.N))
    for ikx in range(self.N):
      kx = ikx*dkx+0.00001 
      for ai,i in zip(filtx,range(1,len(filtx)+1)):    
        xx[ikx,:] += ai*(2*np.cos(i*kx*self.dx)-2)
    return xx

  def Fz(self,filtz):
    nkz = self.N
    dkz = np.pi/(self.dz*(self.N-1))
    zz = np.zeros((self.N,self.N))
    for ikz in range(self.N):
      kz = ikz*dkz+0.00001 
      for ai,i in zip(filtz,range(1,len(filtz)+1)):    
        zz[:,ikz] += ai*(2*np.cos(i*kz*self.dz)-2)
    return zz

  def absk(self):
    kk =  np.zeros((self.N,self.N))
    nkz = self.N
    dkz = np.pi/(self.dz*(self.N-1))
    nkx = self.N
    dkx = np.pi/(self.dx*(self.N-1))
    for ikx in range(self.N):
      kx = ikx*dkx +0.00001 
      for ikz in range(self.N):
        kz = ikz*dkz +0.00001 
        kk[ikx,ikz] = np.sqrt(kx*kx+kz*kz)
    return kk


  def gradient(self,filtx,filtz):
    vp = self.vphase(filtx,filtz)
    dkx = np.pi/(self.dx*(self.N-1))
    dkz = np.pi/(self.dz*(self.N-1))

    v = self.v 
    dt = self.dt
    dx = self.dx
    dz = self.dz 
    kk = self.absk()

    one = np.ones((self.N,self.N))

    DjDvp = self.w*2*(vp-one*v)

    c1 = v*v*dt*dt/(2.0*dx*dx)
    c2 = v*v*dt*dt/(2.0*dz*dz)
    arg = c1*self.Fx(filtx)+\
          c2*self.Fz(filtz)+\
          np.ones((self.N,self.N))

    fix = one/(dt*kk) *(-one/np.sqrt(1-arg*arg))

    gai = [0]*len(filtx)
    gbi = [0]*len(filtz)

    for ai,i in zip(filtx,range(1,len(filtx)+1)):  
      xx = np.zeros((self.N,self.N))
      for ikx in range(self.N):
        kx = ikx*dkx +0.00001         
        xx[ikx,:] = c1*(2*np.cos(i*kx*dx)-2)
      gai[i-1]= (DjDvp*fix*xx).sum()

    for bi,i in zip(filtz,range(1,len(filtz)+1)):  
      zz = np.zeros((self.N,self.N))
      for ikz in range(self.N):
        kz = ikz*dkz +0.00001         
        zz[:,ikz] = c2*(2*np.cos(i*kz*dz)-2)
      gbi[i-1]= (DjDvp*fix*zz).sum()

    return gai,gbi


  def gradientp(self,filtx,filtz):
    from copy import copy
    da = np.abs(np.array(filtx).min())/100000
    ga = [0]*len(filtx)

    f = self.of(filtx,filtz)

    for i in range(len(filtx)):
      ai =  copy(filtx)
      ai[i]+=da
      fai = self.of(ai,filtz)
      ga[i] = (fai-f)/da      

    db = np.abs(np.array(filtz).min())/100000
    gb = [0]*len(filtz)

    for i in range(len(filtz)):
      bi =  copy(filtz)
      bi[i]+=db
      fbi = self.of(filtx,bi)
      gb[i] = (fbi-f)/db 
    return ga,gb




v = 3.
test = LSfilt(0.002,0.01,0.01,v,7,7)
test.setW(0.0)
fx0=[0.1,0,0,0,0]
fz0=[0.1,0,0,0,0]
vphase=test.vphase(fx0,fz0)-v*np.ones((101,101))

fig = plt.figure(1)
ax  = fig.add_subplot(111)
ax.set_xlabel('kx')
ax.set_ylabel('kz')
im = ax.imshow(vphase.T)
fig.colorbar(im)

fx = np.array(fx0)
fz = np.array(fz0)
alpha=1e-9
ofs = []
for it in range(1):
  of = test.of(fx,fz)
  gx,gz = test.gradient(fx,fz)
  gx = np.array(gx)
  gz = np.array(gz)
  
  fx = fx -alpha*gx
  fz = fz -alpha*gz
  if it%10==0:
    print 'it=%d of=%g, |g|=%g '%(it, of,(gx*gx).sum()),fx,fz
  ofs.append(of)

vphase=test.vphase(fx,fz)-v*np.ones((101,101))
# vphase
fig = plt.figure(2)
ax  = fig.add_subplot(111)
ax.set_xlabel('kx')
ax.set_ylabel('kz')
im = ax.imshow(vphase.T)
fig.colorbar(im)

fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.plot(ofs)

plt.show()

