#!/usr/bin/env python

import numpy as np
import rsf.api as rsf
import matplotlib.pyplot as plt


class PointBrowser:
    """
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and previous points
    """
    def __init__(self):
        self.lastind = 0

        self.selected,  = ax.plot([xs[0]], [ys[0]], 'o', ms=12, alpha=0.4,
                                  color='yellow', visible=False)

    def onpress(self, event):
        if self.lastind is None: return
        if event.key not in ('n', 'p'): return
        if event.key=='n': inc = 1
        else:  inc = -1


        self.lastind += inc
        self.lastind = np.clip(self.lastind, 0, len(xs)-1)
        self.update()

    def onpick(self, event):

       if event.artist!=line: return True

       N = len(event.ind)
       if not N: return True

       # the click locations
       x = event.mouseevent.xdata
       y = event.mouseevent.ydata


       distances = np.hypot(x-xs[event.ind], y-ys[event.ind])
       indmin = distances.argmin()
       dataind = event.ind[indmin]

       self.lastind = dataind
       self.update()

    def update(self):
        if self.lastind is None: return

        dataind = self.lastind

        ax2.cla()

        self.selected.set_visible(True)
        self.selected.set_data(xs[dataind], ys[dataind])
        print "using cip %d"%dataind
        print cips.shape
        cip = cips[dataind,:,:,0] 
  
        rht = (nht-1)*dht
        rhx = (nhx-1)*dhx

        ax2.imshow(cip.T,
                   extent=[oht,oht+(nht-1)*dht,ohx,ohx+(nhx-1)*dhx],
                   cmap = plt.get_cmap('gray'),
                   vmin=cips.min()*clipcip,
                   vmax=cips.max()*clipcip,aspect=rht/rhx)
        ax2.set_ylabel(r'$\lambda_x(km)$')
        ax2.set_xlabel(r'$\tau(s)$')
  

        fig.canvas.draw()





def get_axis(File,axis):
  o = File.float("o%d"%axis)
  d = File.float("d%d"%axis)
  n = File.int("n%d"%axis)    
  return o,d,n

def put_axis(File,axis,o,d,n):
  File.put("o%d"%axis,o)
  File.put("d%d"%axis,d)
  File.put("n%d"%axis,n)  





par=rsf.Par()
clip = par.float("pclip",100)/100.0
clipcip = par.float("pclipcip",clip*100)/100.0


Fimage = rsf.Input()
Fpicks = rsf.Input("picks")
Fcip = rsf.Input("cip")

oz,dz,nz = get_axis(Fimage,1)
ox,dx,nx = get_axis(Fimage,2)

image = np.zeros((nx,nz),'f')
Fimage.read(image)

npoints = Fpicks.int("n2")
picks = np.zeros((npoints,2),'f')
print picks.shape
Fpicks.read(picks)

xs = picks[:,0] 
ys = picks[:,1]


X = np.random.rand(100, 200)

fig, (ax, ax2) = plt.subplots(1, 2)
ax.set_title('click on point to plot time series')
ax.imshow(image.T,
          extent=[ox,ox+(nx-1)*dx,oz+(nz-1)*dz,oz],
          cmap = plt.get_cmap('gray'),
          vmin=image.min()*clip,
          vmax=image.max()*clip)
ax.set_xlabel(r'$x(km)$')
ax.set_ylabel(r'$z(km)$')
    


fig2, (ax3) = plt.subplots(1, 1)
ax3.imshow(image.T,
          extent=[ox,ox+(nx-1)*dx,oz+(nz-1)*dz,oz],
          cmap = plt.get_cmap('gray'),
          vmin=image.min()*clip,
          vmax=image.max()*clip)
ax3.set_xlabel(r'$x(km)$')
ax3.set_ylabel(r'$z(km)$')

 
line, = ax.plot(xs, ys, '.', picker=5)  # 5 points tolerance




ohz,dhz,nhz = get_axis(Fcip,1)
ohx,dhx,nhx = get_axis(Fcip,2)
oht,dht,nht = get_axis(Fcip,3)

print nhz,nhx,nht
cips = np.zeros((npoints,nht,nhx,nhz),'f')

Fcip.read(cips)






#
browser = PointBrowser()
#
fig.canvas.mpl_connect('pick_event', browser.onpick)
fig.canvas.mpl_connect('key_press_event', browser.onpress)




plt.show()

Fimage.close()
Fpicks.close()
Fcip.close()


