import numpy as np
import matplotlib
matplotlib.rcParams.update({'font.size':14})
import matplotlib.pyplot as plt

def sin(theta):
  return np.sin(theta*np.pi/180.)

def cos(theta):
  return np.cos(theta*np.pi/180.)

def tan(theta):
  return np.tan(theta*np.pi/180.)

def arctan(x):
  return (180./np.pi)*np.arctan(x)



class TI_medium:
  def __init__(self,vp0,vs0,eps,delta,mode='fast'):
    s = {'fast':1.,'slow':-1}

    self.vp0 = vp0
    self.vs0 = vs0
    self.eps = eps
    self.delta = delta
    self.f = 1. - np.power(vs0/vp0,2.)
    self.mode = s[mode]

  def phase_velocity(self,theta):
    a = 1. + self.eps*np.power(sin(theta),2) -0.5*self.f
    b = self.mode*self.f*0.5
    c = 1. + 2.*self.eps*(1./self.f)*np.power(sin(theta),2.)
    c *=c 
    d = 2*(self.eps-self.delta)*np.power(sin(2*theta),2.)/self.f
    b *= np.sqrt(c-d)
    return self.vp0*np.sqrt(a +b)

  def dphase_velocity(self,theta,h=0.001):
    f = self.f
    eps = self.eps  
    delta = self.delta
    vp0 = self.vp0
    vs0 = self.vs0 
  
    vp = self.phase_velocity(theta)
    sq = np.sqrt((1. + (4./f)*(sin(theta)*sin(theta))*(2*delta*cos(theta)*cos(theta) - eps*cos(2.*theta))
           + (np.power((2*eps/f),2))*(np.power(sin(theta),4))))

    t1 = (np.power(vp0,2.)*0.5*eps*sin(2*theta))/vp
    t2 = (0.5*sin(2.*theta)*(2*delta*(cos(theta)*cos(theta)) - eps*cos(2.*theta)) + sin(theta)*sin(theta)*(eps*sin(2.*theta) 
          - delta* sin(2.*theta)) + (2/f)*np.power(eps,2.)*(np.power(sin(theta),3.)*cos(theta)))
    return t1 +self.mode*((vp0*vp0)/vp)*(t2/sq)

  def phase2group(self,theta):
    t = tan(theta)
    v = self.phase_velocity(theta)
    dv = self.dphase_velocity(theta)
    ivdv = (1./v)*dv
    tantsi = (t + ivdv)/(1.-t*ivdv)
    return arctan(tantsi)

  def vphase2vgroup(self,theta):
    v = self.phase_velocity(theta)
    dv = self.dphase_velocity(theta)
    ivdv = (1./v)*dv
    return v*np.sqrt(1+(ivdv*ivdv))    
    


class layer:
  def __init__(self,thickness,vp0,vs0,eps,delta,mode='fast'):
    self.t = TI_medium(vp0,vs0,eps,delta,mode)
    self.d = thickness


    


def ray(layers,theta=20.0):
  t = 0.0 

  nl = len(layers)
  rpath = np.zeros((2*nl+1,2))

  layer = layers[0]
  layer2 = layers[1]
  tsi = layer.t.phase2group(theta)
  m = tan(tsi)
  y1 = layer.d +rpath[0,1]
  x1 = m*layer.d +rpath[0,0]
  rpath[1,0] = x1
  rpath[1,1] = y1

  t += (layer.d/cos(tsi))/layer.t.vphase2vgroup(theta) # time first layer downward

  theta2 = snell(theta,layer,layer2)
  tsi = layer2.t.phase2group(theta2)
  t += (layer2.d/cos(tsi))/layer2.t.vphase2vgroup(theta2) # time second layer downward

  y2= layer.d+layer2.d 
  m = tan(tsi)
  x2 = m*layer2.d +rpath[1,0] 
  rpath[2,0] = x2
  rpath[2,1] = y2

  # reached bottom layer
  x3 = m*layer2.d + rpath[2,0]
  y3 = layer.d
  rpath[3,0] = x3
  rpath[3,1] = y3

  theta = snell(theta2,layer2,layer)
  tsi = layer.t.phase2group(theta)
  m = tan(tsi)
  x4 = m*layer.d +rpath[3,0]
  y4 = 0.0

  rpath[4,0] = x4
  rpath[4,1] = y4

  return (rpath,2.*t)

def snell(theta1,layer1,layer2):
  from scipy.optimize import fsolve,fminbound
  v1 = layer1.t.phase_velocity(theta1)
  sx1 = sin(theta1)/v1
  theta2 = fminbound(f,-90.,90.,(layer2,sx1))
  return theta2

def f(x,layer,sx):
  m = sx - sin(x)/layer.t.phase_velocity(x) 
  return np.abs(m)


