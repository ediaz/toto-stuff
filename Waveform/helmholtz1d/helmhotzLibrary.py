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

  def backwardDiff(self):
    md = np.ones(self.n)
    md[0] = 0
    ld = np.ones(self.n)
    return sparse.spdiags([md,-ld],[0,-1],self.n,self.n)

  def forwardDiff(self):
    return -1*self.backwardDiff().T

   


hm = helmhotz1D(0,1,4)

bd = hm.backwardDiff().todense()
fd = hm.forwardDiff().todense()

print bd
print fd

print fd*bd
