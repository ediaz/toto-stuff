import numpy as np
import scipy.sparse as sparse
from scipy.sparse import linalg as splinalg
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

class helmhotz1D:
  '''
  This class solves:
     1               d  1   d
   ------- w^2 u  + -- ---  -- u = -s
   rho v^2          dx rho  dx
  '''
  def __init__(self,o,d,n):
    self.o = o
    self.d = d
    self.n = n

  def setVel(self,vel):
    if len(vel) != self.n:
      print "error in dimensions"
    else: 
      self.vel = vel

  def setRho(self,rho):
    if len(vel) != self.n:
      print "error in dimensions"
    else: 
      self.rho = rho*(1.0+.0j)

  def setSource(self,sou):
    self.f = 1.0*sou

  def d2dt(self,w):
    k2 = 1.*w*w/(self.rho*self.vel*self.vel)
    return sparse.spdiags([k2],[0],self.n,self.n)    
    
  def dxRhoDx(self):
    irho = sparse.spdiags([1./self.rho],[0],self.n,self.n)
    return self.forwardDiff()*irho*self.backwardDiff()

  def lhs(self,w):
    return self.d2dt(w)+self.cden()  

  def cden(self):
    ud = np.ones(self.n)
    md = -2*np.ones(self.n)
    ld = np.ones(self.n)
    return sparse.spdiags([ud,md,ld],[1,0,-1],self.n,self.n)*1./self.d**2

  def backwardDiff(self):
    md = np.ones(self.n)
    ld = np.ones(self.n)
    return sparse.spdiags([md,-ld],[0,-1],self.n,self.n)*1./self.d

  def forwardDiff(self):
    return -1*self.backwardDiff().T

  def solve(self,w):
    A = self.lhs(w)
    b = self.f
    return splinalg.spsolve(A,b)
################################################################### 

def source(o,d,n,x,w):
  s = np.zeros(n,dtype=complex)
  s[int((x-o)/d)] = -1.0+0.j
  return s


def complexvel(vel,Q):
  v = vel - 1j*vel/(2.*Q)
  return v


o,d,n = 0.0,0.01,101
x = d 

s = source(0,d,n,x,1.)
sign = (np.sign(np.arange(-n/2,n/2,1)+0.0001)+1)*0.5
varQ = (1-sign)*10000+sign*10
vel = np.ones(n)*1.5 #3*(np.sign(np.arange(-n/2,n/2,1)+0.0001)*0.5+2)

velQ = complexvel(vel,varQ)

rho = np.ones(n) #vel/1.5#*(np.sign(np.arange(-n/2,n/2,1)+0.0001)*0.5+2) 

hm = helmhotz1D(o,d,n)
hm.setVel(vel)
hm.setRho(rho)
hm.setSource(s)

nt = 100
dt = 1./30
df = 1./(dt*nt) # 1/Total time 

print 'dt=',dt,'df=',df



twfl = np.zeros((nt,n),dtype=complex)
tau = 1. 
iw=0
for f in np.arange(df,nt*df,df):
  iw +=1
  w = f*2*np.pi
  wfl = hm.solve(w+1.j/tau)
  for it in range(nt):
    twfl[it,:] += wfl*np.exp(-iw*iw/(2.*(30*30)))*np.exp(-w*1j*it/nt)

twfl = twfl.real/nt
for it in range(nt):
  tt = 0.35*it*dt/tau
  twfl[it,:] = twfl[it,:]*np.exp(tt)

xv = np.arange(o,n*d+o,d)
fig = plt.figure(111)
for t in range(10,31,2):
  ax = fig.add_subplot(111)
  ax.plot(xv,twfl[t,:])
  plt.ylim(ymax = twfl.max(), ymin = -twfl.max())
#

#wfl1 = hm.solve(2*np.pi*(20))
#hm = helmhotz1D(o,d,n)
#hm.setVel((1-sign)*1.5+sign*2)
#hm.setRho(rho)
#hm.setSource(s)
#wfl2 = hm.solve(2*np.pi*(20))
#

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(xv,wfl.real)
ax1.set_xlabel("x(km)")
ax2 = fig.add_subplot(212)
ax2.plot(xv,wfl.real)
ax2.set_xlabel("x(km)")


##
plt.show()
