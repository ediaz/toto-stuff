import numpy as np
import scipy.sparse as sparse
from scipy.sparse import linalg as splinalg


class helmhotz1D:
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
      self.rho = rho

  def setSource(self,sou):
    self.f = -1*sou

  def d2dt(self,w):
    k2 = 1.*w*w/(self.rho*self.vel*self.vel)
    return sparse.spdiags([k2],[0],self.n,self.n)    
    
  def dxRhoDx(self):
    irho = sparse.spdiags([1./self.rho],[0],self.n,self.n)
    return self.forwardDiff()*irho*self.backwardDiff()

  def lhs(self):
    return self.d2dt()+dxRhoDx()  



  def backwardDiff(self):
    md = np.ones(self.n)
    md[0] = 0
    ld = np.ones(self.n)
    return sparse.spdiags([md,-ld],[0,-1],self.n,self.n)*1./self.d

  def forwardDiff(self):
    return -1*self.backwardDiff().T

   


hm = helmhotz1D(0,1,4)

bd = hm.backwardDiff().todense()
fd = hm.forwardDiff().todense()



print "backward der"
print bd
print "forward der"
print fd
print "second derivative"
print fd*bd
