import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import scipy.linalg as linalg
import matplotlib
matplotlib.rcParams.update({'font.size': 22})



def getalpha(dmod,fmap,feval,df,m,fp,pk,maxiter):
  '''
  dmod: forward modeled data at current model
  fmap: forward mapping function (takes model m as input)
  feval: given data, computes the of value
  df : pointer to gradient function
  m  : current model
  fp : current ASCENT direction
  maxiter : how many max functions evaluation you
            want to use to find opt alpha? 
            normally only 1 is used 
  '''
  a = 0.1
  dma = fmap(m-a*fp)
  r = (dmod-dma)/a
  v = np.dot(r,r)

  alpha = 0.5*np.dot(pk,pk)/v # alpha equation: page 62, chapter 4 from his notes
                              # magick factor 0.5?
  i=0  
  return alpha,i


def steepest_descent(fmap,feval,df,m,tol=1.0e-15,maxiter=30,optalpha=6,ax=None):
  '''
  fmap : pointer to function that computes data
  feval:  pointer to the objective function that given data computes the of
  df   :  pointer to a function that computes the gradient of 
          feval at model m
  m    :  initial model
  tol  :  tolerance, the solver quits if the of did not decrease
          more than tol
  maxiter : maximum number of iterations
  ax   : optionally input for matplotlib axes
  optalpha: maximum number of extra evaluations to search for optimum alpha
  '''
  fo = 0
  finitial = feval(m)
  for it in range(maxiter):
    dmod = fmap(m) 
    fm = feval(dmod)
    fp = df(m)

    #alpha, new_fval, new_x,i = simple_line_search_owlqn(feval, fm, m, -fp)
    alpha,i = getalpha(dmod,fmap,feval,df,m,fp,fp,optalpha)
    if ax: 
      ax.quiver(m[0],m[1],-alpha*fp[0],-alpha*fp[1],scale_units='xy', scale=1,color='k')
    m -= alpha*fp
    print 'of= %15.8g'%fm,'feval= %3d'%(i+2),'step=',alpha
  
    if abs(fm-fo) < tol or fm <2*tol:
      break
    fo = fm
  return m


def NLconjugate_gradient(fmap,feval,df,m,tol=1.0e-15,maxiter=30,optalpha=6,ax=None):
  '''
  fmap : pointer to function that computes data
  feval:  pointer to the objective function that given data computes the of
  df   :  pointer to a function that computes the gradient of 
          feval at model m
  m    :  initial model
  tol  :  tolerance, the solver quits if the of did not decrease
          more than tol
  maxiter : maximum number of iterations
  ax   : optionally input for matplotlib axes
  optalpha: maximum number of extra evaluations to search for optimum alpha
  '''
  dmod = fmap(m)
  fo = feval(dmod)
  fp = df(m) 
  alpha,i = getalpha(dmod,fmap,feval,df,m,fp,fp,optalpha)
  if ax:
    ax.quiver(m[0],m[1],-alpha*fp[0],-alpha*fp[1],scale_units='xy', scale=1,color='g')
  print 'of= %15.8g'%fo,'feval= %3d'%(i+2),'step=',alpha
  m -= alpha*fp
  so = -fp
  dxo = so
  fm = 0. 
  early = False
  for it in range(maxiter-1):
    dmod = fmap(m)
    fm = feval(dmod)
    dxn = -df(m)
    beta = np.dot(dxn,dxn)/np.dot(dxo,dxo) #Fletcher-Reeves
    sn = dxn+beta*so
    alpha,i = getalpha(dmod,fmap,feval,df,m,-sn,-dxn,optalpha)
    so = sn
    if ax:
      ax.quiver(m[0],m[1],alpha*sn[0],alpha*sn[1],scale_units='xy', scale=1,color='g')
    print 'of= %15.8g'%fm,'feval= %3d'%(i+2),'step=',alpha
    m += alpha*sn
    dxo = dxn  
    if abs(fm-fo) < tol:
      early=True
      break
    fo = fm
  if early: print "finished early! only %d iterations were used"%(it+1)
  return m










def fmap((m1,m2)):
  '''
  forward map function: 
    given a model m produces data dmod
  '''
  return np.array((m1+m2,-2.*m1+3.*m2))


def hessian((m1,m2)):
  a11 = 10
  a12 = -10
  a21 = -10
  a22 = 20
  return np.matrix([[a11,a21],[a12,a22]])





def data((m1,m2)):
  '''
  this function returns the observed data
  it looks fancy because I had to adapt it
  to acept vectorized operations
  '''
  a = np.ones(m1.shape)
  return np.array((2.*a,1.*a))

def of(dmod):
  '''
  give modeled data computes the objective function
  '''
  dobs = data(dmod)
  res = dmod - dobs
  return np.sum(res*res,axis=0)

def grad((m1,m2)):
  '''
  returns the gradient vector at the 
  model m=(m1,m2) 
  '''
  dm1 = 2.*(m1+m2-2.)-4.*(3.*m2-2.*m1-1)
  dm2 = 2.*(m1+m2-2.) +6.*(3.*m2-2.*m1-1)
  g = np.array((dm1,dm2))
  return g



mo= np.array((6.,5.))
m = np.zeros(2)
m = mo

A = np.array([[1,1.],[-2,3]])
b = np.array([2,1.])


# here I produce the contour plot of the of:
x = y = np.arange(-0.,7,0.1)
X,Y = np.meshgrid(x,y)
dmod = fmap((X,Y))
Z = of(dmod)

fig = plt.figure()
ax = fig.add_subplot(111,aspect='equal')
ax.contour(X,Y,Z,60,colors='k')
ax.plot(1.,1.,'yo',markersize=10)
ax.set_xlabel("$m_1$")
ax.set_ylabel("$m_2$")
ax.set_title("toy problem: figure 17")




m= np.array((6.,5.))
ax.plot(m[0],m[1],'bo',markersize=10)






mu = steepest_descent(fmap,of,grad,m,tol=1.0e-20,maxiter=31,optalpha=6,ax=ax)
#print mu

m= np.array((6.,5.))
mu = NLconjugate_gradient(fmap,of,grad,m,tol=1.0e-20,maxiter=4,optalpha=6,ax=ax)
#print mu

plt.show()
