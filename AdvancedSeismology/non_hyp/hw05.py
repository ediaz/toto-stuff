import numpy as np
import matplotlib
matplotlib.rcParams.update({'font.size':14})
import matplotlib.pyplot as plt




xt  = np.loadtxt('PwaveNonHyp')
x = xt[:,0]
t = xt[:,1]

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(x,t)
ax.invert_yaxis()
ax.set_xlabel("x(km)")
ax.set_ylabel("t(s)")
ax.set_title('observed traveltimes')



class layer:
  def __init__(self,tick,Vp0,eps,delta):
    self.tick = tick
    self.vp0 = Vp0
    self.eps = eps
    self.delta = delta
    self.eta = (eps-delta)/(1.+(2.*delta))
    self.t0 = 2*tick/Vp0
    self.vnmo = self.vp0*np.sqrt(1.+2.*delta)
    self.vhor = self.vnmo*np.sqrt(1+2*self.eta)


layer1 = layer(1.480,1.941,0.230,0.011)
layer2 = layer(1.980-1.480,2.168,0.247,0.010)
layer3 = layer(2.244-1.980,2.685,0.264,0.016)


# effective parameters:
t0 = 0. 
Vnmo = 0. 
eta = 0. 
for layers in [layer1,layer2,layer3]:
  t0 += layers.t0
  Vnmo += np.power(layers.vnmo,2.)*layers.t0
  eta += np.power(layers.vnmo,4.)*(1.+8.*layers.eta)*layers.t0


Vnmo = np.sqrt(Vnmo/t0)
eta = 1./8.*(1./(np.power(Vnmo,4.)*t0)*eta-1.)
Vhor = Vnmo*np.sqrt(1.+2.*eta)
print t0,Vnmo,eta,Vhor


def tquartic(x,Vnmo2=Vnmo,Vhor2=Vhor,C=1.):
  p1 = t0*t0 + x*x/(Vnmo2*Vnmo2)
  p2 = (Vhor2*Vhor2-Vnmo2*Vnmo2)*np.power(x,4.)
  p3 = ((t0*t0*np.power(Vnmo2,4.)) + (C*np.power(Vhor2,2.)*(x*x)))*(Vnmo2*Vnmo2);
  t2  = np.sqrt(p1 - p2/p3)
  return t2

N=1001
l2 = np.zeros(N)
C = np.linspace(0.9,1.4,N)
for i in range(N):
  t2 = tquartic(x,C=C[i])
  l2[i] = np.power((t-t2),2.).sum() 

fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.plot(C,l2)
ax.set_xlabel("C")
ax.set_ylabel("L2 error")




# optimum C:
optC = C[np.argmin(l2)]
print "optimum C=", optC,l2[np.argmin(l2)]

fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.plot(x,t)
ax.plot(x,tquartic(x,C=optC))
ax.set_xlabel("x(km)")
ax.set_ylabel("Time(s)")
ax.invert_yaxis()




nvnmos = 201
nvhors = 201
vnmos = np.linspace(1.7,2.5,nvnmos)  #Vnmo*0.7,Vnmo*1.3,nvnmos)
vhors = np.linspace(2.0,2.8,nvhors)  #Vhor*0.7,Vhor*1.2,nvhors)

l2vel = np.zeros((nvnmos,nvhors))

for ivnmo in range(nvnmos):
  vnmo = vnmos[ivnmo]
  for ivhor in range(nvhors):
    vhor = vhors[ivhor]
    t2 = tquartic(x,vnmo,vhor,C=1.)
    err = np.power((t2-t),2.).sum() 
    l2vel[ivnmo,ivhor] = err


(optivnmo,optivhor) =  np.unravel_index(l2vel.argmin(), l2vel.shape)

optvnmo = vnmos[optivnmo]
optvhor = vhors[optivhor]



print optvnmo,optvhor

effeta = (np.power(optvhor/optvnmo,2.)-1.)*0.5
print "searched effective eta:",effeta,eta

t3 = tquartic(x,optvnmo,optvhor,C=1.)
ax.plot(x,t3)

fig = plt.figure(4)

X,Y = np.meshgrid(vnmos,vhors)

ax = fig.add_subplot(111)
cs = ax.contour(X,Y,l2vel.T,180)
ax.clabel(cs, inline=1, fontsize=10)

ax.set_xlabel("Vnmo(km/s)")
ax.set_ylabel("Vhor(km/s)")



fig = plt.figure(5)
ax = fig.add_subplot(111)
ax.plot(x,t3-t,'r')
ax.plot(x,tquartic(x,C=optC)-t)
ax.set_xlabel("x(km)")
ax.set_ylabel("Error(s)")



plt.show()
