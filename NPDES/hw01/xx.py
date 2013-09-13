import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 22,'legend.fontsize': 22})




def theta_w(a,s,theta):
    if(theta >a):
      w2 = 1.0;
    else:
      w2 = np.exp(-((theta-a)*(theta-a))/(2.0*s*s));
    return w2;


a = 20.0
b = 70.0

PI = np.arccos(-1.0)


a = np.cos(a*PI/180.0) 
b = np.cos(b*PI/180.0) 
s = np.abs(a-b)/3.0


x = np.zeros(180)
x2 = np.zeros(180)
y = np.zeros(180)
for it in range(180):
  x[it] = np.cos(it*PI/180)
  x2[it] = (it)
  y[it] = theta_w(a,s,x[it])

fig = plt.figure(1)
ax1 = plt.subplot(111)

ax1.plot(x2,y,linewidth=5)
ax1.set_xlabel(r"$\theta(^\circ)$")
ax1.set_ylabel(r"$W(\theta)$")

fig = plt.figure(2,figsize=(8,3))
ax1 = plt.subplot(111)

ax1.plot(x,y,linewidth=5)
ax1.set_xlabel(r"$\cos(\theta)$")
ax1.set_ylabel(r"$W(\theta)$")

plt.show()
 
