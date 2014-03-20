import numpy as np
import scipy.sparse as sparse
from scipy.sparse import linalg as splinalg
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import scipy.linalg as linalg
import matplotlib
plt.rcParams.update({'font.size': 20,'legend.fontsize': 20})



def getalpha(fo,feval,df,m,fp,maxiter):
  '''
  fo : of value at current model
  feval : pointer to of
  df : pointer to gradient function
  m  : current model
  fp : current ascent direction
  maxiter : how many max functions evaluation you
            want to use to find opt alpha? 
  '''
  a = 0.02
  fpa = df(m-a*fp)
  J = -(fpa-fp)/a 
  n = (fp*fp).sum(); d = (J*fp).sum()+1.e-30  
  alpha =(n/d)
  dm = -alpha*fp
  f = 1.0 
  for i in range(maxiter):
    if feval(m +f*dm)<fo:
      break
    else:
      f *=0.5
  return alpha*f,i


def steepest_descent(feval,df,m,tol=1.0e-15,maxiter=30,optalpha=6,ax=None):
  '''
  feval:  pointer to the objective function that takes m 
  df   :  pointer to a function that computes the gradient of 
          feval at point m
  m    :  initial model
  a    :  base step size to compute f(m+a*df)
  tol  :  tolerance, the solver quits if the of did not decrease
          more than tol
  maxiter : maximum number of iterations
  ax   : optionally input for matplotlib axes
  optalpha: maximum number of extra evaluations to search for optimum alpha
  '''
  fo = 0
  finitial = feval(m)
  for it in range(maxiter):
  
    fm = feval(m)
    fp = df(m)

    #alpha, new_fval, new_x,i = simple_line_search_owlqn(feval, fm, m, -fp)
    alpha,i = getalpha(fm,feval,df,m,fp,optalpha)
    if ax: 
      ax.quiver(m[0],m[1],-alpha*fp[0],-alpha*fp[1],scale_units='xy', scale=1)
    m -= alpha*fp
    print 'of= %15.8g'%fm,'feval= %3d'%(i+2),'step=',alpha
  
    if abs(fm-fo) < tol:
      break
    fo = fm
  return m




def of((m1,m2)):
  return 4*m1*m1 -4*m1*m2 +2*m2*m2 

def grad((m1,m2)):
  dm1 = 2*np.sin(m1)*np.sin(m2)*np.cos(m1)*np.sin(m2) 
  dm2 = 2*np.sin(m1)*np.sin(m2)*np.cos(m2)*np.sin(m1)
  return np.array((dm1,dm2))



x = y = np.arange(-2.,2,0.05)
X,Y = np.meshgrid(x,y)
zs = np.array([of((x,y)) for x,y in zip(np.ravel(X), np.ravel(Y))])
Z = zs.reshape(X.shape)


fig = plt.figure()
ax = fig.add_subplot(111,aspect='equal')
cs = ax.contour(X,Y,Z,60)
ax.clabel(cs, inline=1, fontsize=10)



mo = (1.,-.8)
ax.plot(mo[0],mo[1],'yo')
ax.set_xlabel("$m_1$")
ax.set_ylabel("$m_2$")
#
#
#
#x = y = np.arange(0.,7.01,0.5)
#X,Y = np.meshgrid(x,y)
#xx = -grad((X,Y))
#U = xx[0]
#V = xx[1]
##ax.quiver(X,Y,U,V)
#
#
#m= np.array(mo)
#mu = steepest_descent(of,grad,mo,tol=1.0e-11,maxiter=30,optalpha=6,ax=ax)
##print mu
#
##NLconjugate_gradient(of,grad,m,ax=ax)
#

plt.show()
