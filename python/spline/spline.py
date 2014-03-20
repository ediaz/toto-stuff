import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

#####################################################


nx,nz = 231,76
n = nx*nz

x = np.zeros(n,'f')
y = np.zeros(n,'f')
z = np.zeros(n,'f')


n2 = len(range(0,nx,4))*len(range(0,nz,4))
xx = np.zeros(n2,'f')
yy = np.zeros(n2,'f')




marm = np.fromfile('marmousi.dat','f')
marm = np.reshape(marm,(nx,nz))


i = 0 
for i2 in range(nx):
  x[i] = i2*1.
  for i1 in range(nz):
    y[i] = i1*1.
    z[i] = marm[i2,i1]
    i += 1
i=0
for i2 in range(0,nx,4):
  xx[i] = i2
  for i1 in range(0,nz,4):
    yy[i] = i1
    i+=1

[tck] = interp.bisplrep(x, y, z, kx=3, ky=3, eps=1e-16,full_output=1, nxest=None, nyest=None, quiet=1)


fig = plt.figure(1,figsize=(10,5))
ax  = fig.add_subplot(111)







ax.imshow(marm.T)

plt.show()
