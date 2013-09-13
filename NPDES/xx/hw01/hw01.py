import numpy as np
import heatEquation as heq
from functions_hw01 import *
from plotting_hw01 import *


  
print "%2s&%10s&%14s&%14s\\\\"%('dx','dt','error','EOC')

ntest = 5

xh = np.zeros(ntest)
eh = np.zeros(ntest)
ifig=1
for k in range(ntest):
  dx = dt = 0.1/2**k
  ox = 2.0
  nx = int((3.0+1e-10-ox)/dx)+1
  ot = 0.0
  nt = int((2.0)/dt)+1
  sx = (ox,dx,nx)
  st = (ot,dt,nt)
  
  ue =  uexact(sx,st)
  a = heq.CrankNicolson(f,v,g1,g2,st,sx)
  sol = a.solve()
  e = compute_error(sol,sx,st)
  ue = uexact(sx,st)
  diff = sol-ue

  display2d_3(nt,ot,dt,nx,ox,dx,sol,ue,diff,ifig,"",'report/fig/solutions_k%d.png'%k)
  xh[k] = dx
  eh[k] = e  
  print '%-10.5g&%-10.5g&%-15.6e&%-15.4e \\\\'%(dx,dt,e,2*dx*dx)
  ifig+=1

error_plot(xh,eh,ifig,'report/fig/error.png')

# show() # uncomment to see the plots. If is commented the plots are directly written to png files
