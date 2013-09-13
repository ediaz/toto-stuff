import numpy as np
import heatEquation as heq
from functions_hw02 import *
from plotting_hw02 import *

##########
N = 101
M = N 
ox,dx,nx = 0,1./(N-1),N
ot,dt,nt = 0,1./(M-1),M

##########

sx = (ox,dx,nx)
st = (ot,dt,nt)

obj = heq.t3level(f1,v1,st,sx,a=1./2.)
u3l = obj.solve()

obj = heq.implicit(f1,v1,st,sx,a=1./2.)
ube = obj.solve()
#Figure 1 in the report
display2d_2(nt,ot,dt,nx,ox,dx,u3l,ube,fn=1,title="",save=None)

obj = heq.t3level(f2,v2,st,sx,a=1./12.)
u3l = obj.solve()

ue = uexact(sx,st)
#Figure 2 in the report
display2d_3(nt,ot,dt,nx,ox,dx,u3l,ue,ue-u3l,fn=4,title="",save=None)

test = [50,41,81,161,321]
xh = np.zeros(len(test))
eh = np.zeros(len(test))

print "%2s%10s%14s        N"%('dx','dt','error')
for i in range(len(test)):
  M = N =  test[i] 
  ox,dx,nx = 0,1./(N-1),N
  ot,dt,nt = 0,1./(M-1),M
  ##########

  sx = (ox,dx,nx)
  st = (ot,dt,nt)

  
  obj = heq.t3level(f2,v2,st,sx,a=1./12.)
  u3l = obj.solve()
  ue = uexact(sx,st)

  
  l2 =  l2err(ue-u3l)
  print '%-10.5g%-10.5g%-15.6e%2d'%(dx,dt,l2,N)
  xh[i] = dx
  eh[i] = l2 

error_plot(xh,eh,5)
print np.polyfit(np.log(xh),np.log(eh),1)
show()
