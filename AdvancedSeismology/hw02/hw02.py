import numpy as np
import scipy.sparse as sparse
from scipy.sparse import linalg as splinalg
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm



def VpApprox(vp0,eps,delta,vp0vs0,theta):
  f = 1.-np.power(1./vp0vs0,2)
  sint = np.sin(theta)
  sin2t = sint*sint
  sin4t = np.power(sint,4)
  cost = np.cos(theta)
  cos2t = cost*cost
  
  v = vp0*np.sqrt(1.+
      2*delta*sin2t*cos2t+
      2*eps*sin4t +
      4./f*(eps*sin2t +delta*cos2t)*sin4t*cos2t
      )

  v2 = vp0*vp0*(
      1. + 2*delta*sin2t*cos2t
      + 2*eps*sin4t +
      4/f*(eps-delta)*(eps*sin2t+delta*cos2t)*sin4t*cos2t
      )
  return np.sqrt(v2)


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
        -0.5*f +
        0.5*f*np.sqrt(a*a-b))

  return np.sqrt(v2)



def goQuestion1():
  thetad = np.arange(0.,90.,1.)
  theta = thetad*np.pi/180
  vp0vs0 = 2.
  eps = 0.6
  deltas = [0.1,0.3,0.5]
  vp0 = 1.
  ifn = 2
  for delta in deltas:
    vpe = VpExact(vp0,eps,delta,vp0vs0,theta)
    vpa = VpApprox(vp0,eps,delta,vp0vs0,theta)
    fig = plt.figure(ifn)
    ax = fig.add_subplot(111)
    ax.plot(thetad,vpe,label='exact phase velocity',linestyle='--')
    ax.plot(thetad,vpa,label='approximated phase velocity')
    ax.legend(loc='upper left')
    ax.set_xlabel("Phase angle($^\circ$)")
    ax.set_ylabel("Phase velocity")
    #plt.savefig('report/Fig/q1_delta%d'%(int(delta*10)),bbox_inches='tight')
    ifn+= 1


def goQuestion2():
  thetad = np.arange(0.,90.,1.)
  theta = thetad*np.pi/180
  vp0 = 1.
  d1 = {'eps':0.6,'delta':0.1}
  d2 = {'eps':0.6,'delta':0.5}
  m = [d1,d2]
  vp0vs0s = np.arange(1.5,4.,2.5/4)

  ifn = 10
  for i in range(2):
    mi = m[i]
    eps = mi['eps']
    delta = mi['delta']
    mtag = 'model%d'%(i+1)
    fig = plt.figure(ifn)
    sym = ['ko','k-','k^','kd']
    fig = plt.figure(ifn,dpi=280)
    ax = fig.add_subplot(111) 
    for vp0vs0 in vp0vs0s:
      vpe = VpExact(vp0,eps,delta,vp0vs0,theta)
      ax.plot(thetad,vpe,label='vp0/vs0 = %3.2g, eps=%3.2g, delta=%3.2g'%(vp0vs0,eps,delta))

    ax.legend(loc='upper left')
    ax.set_xlabel("Phase angle($^\circ$)")
    ax.set_ylabel("Phase velocity")
    plt.savefig('report/Fig/q2_model%d'%(i+1))
    ifn+= 1

    vp0vs0x = np.arange(1.5,4.,0.1)
    thetas = 45*np.pi/180
    vpe = VpExact(vp0,eps,delta,vp0vs0x,thetas)
    fig = plt.figure(20+i)
    ax = fig.add_subplot(111) 
    ax.set_xlabel("Vp0/Vs0(theta=45)")
    ax.plot(vp0vs0x,vpe)
    plt.savefig('report/Fig/dependance%d'%(i+1))
  
  


goQuestion2()
#goQuestion1()
plt.show()
