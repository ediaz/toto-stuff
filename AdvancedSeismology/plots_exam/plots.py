import numpy as np
import scipy.sparse as sparse
from scipy.sparse import linalg as splinalg
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm



def vsv(sigma,vs0=1.4):
  theta = np.arange(0,90,0.1)*np.pi/180.
  v = np.zeros(len(theta))
  for i in range(len(theta)):
    t = theta[i]
    v[i] = vs0*(1+sigma/4.*np.power(np.sin(2*t),2))
  return (v,theta)  


def VpExact(vp0,eps,delta,vp0vs0,theta):
  f = 1.-np.power(1./vp0vs0,2)
  sint = np.sin(theta)
  sint2 = np.sin(2*theta)
  sin2t2 = sint2*sint2
  sin2t = sint*sint
  sin4t = np.power(sint,4)
  cost = np.cos(theta)
  cos2t = cost*cost
  
  a = (1+ 2.*eps*sin2t/f)
  b = 2.*(eps-delta)*sin2t2/f
  v2 = vp0*vp0*(1. 
        + eps*sin2t 
        -0.5*f -
        0.5*f*np.sqrt(a*a-b))

  return np.sqrt(v2)


  
fig = plt.figure(1)
ax  = fig.add_subplot(111)
v,t = vsv(0.48)
ax.plot(t,v)

fig2 = plt.figure(2)
ax2  = fig2.add_subplot(111)

for sigma in np.arange(0.48,1.5,.2):
  v,t = vsv(sigma)
  vv = (1./v)*np.cos(t)
  vh = (1./v)*np.sin(t)
  ax2.plot(vh,vv)
ax2.set_title("sv slowness surface")


plt.show()
