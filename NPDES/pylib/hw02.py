import numpy as np
import heatEquation as heq
from graphRecipes import *
import hw01_functions as fhw01

    

M = 501
N = 11
ot,dt,nt = fhw01.ganesh2Toto(1.,M)
ox,dx,nx = fhw01.ganesh2Toto(1.,N)



fcns = {'f1':fhw01.f1,'f2':fhw01.f2,'f3':fhw01.f3,
        'v1':fhw01.v1,'v2':fhw01.v2,'v3':fhw01.v3,
        'n1_v':"$v(x) = (1-{sign}(x-0.5))x+(1+{sign}(x-0.5))(1-x)$",
        'n2_v':"$v(x) = \sin(2\pi x)$",
        'n3_v':"$v(x) = x^2-x$",
        'n1_f':"$f(x,t) = \sin(\pi x)\exp(-2t)$",
        'n2_f':"$f(x,t) = \sin(\pi x)\exp(t)$",
        'n3_f':"$f(x,t) = 64(x(1-x)t)-4x$"}




################## Excercise 2 ######################
fnum=1
#for p in range(1,4,1):
#  f = fcns['f%d'%p]
#  v = fcns['v%d'%p]
#  fn = fcns['n%d_f'%p]
#  vn = fcns['n%d_v'%p]
#  title = fn+'\n'+vn+'\n'
#
#  print "------------\nsolving explicit Euler with:\n"+title+"------------"
#  a1 = heq.implicitEuler(f,v,\
#                       nt,ot,dt,\
#                       nx,ox,dx,\
#                       a=1.0)
#  sol1 = a1.solve()
#  display2d(nt,ot,dt,nx,ox,dx,sol1,fnum,title)
#  j = int(0.2/dt)
#  graph_A(sol1,nx,ox,dx,j,fnum+1,dt,title)
#  fnum+=2
#

############### Exercise 3 #######################
print "Exercise 2.3: controled solution"
exact = fhw01.u_exact(nt,nx)
display2d(nt,ot,dt,nx,ox,dx,exact,fnum,"Exact solution")
a1 = heq.implicitEuler(fhw01.f_exact,fhw01.v_exact,\
                       nt,ot,dt,\
                       nx,ox,dx,\
                       a=1.0)

sol1 = a1.solve()
display2d(nt,ot,dt,nx,ox,dx,sol1,fnum+1,"eE solution")

e = exact - sol1
display2d(nt,ot,dt,nx,ox,dx,e,fnum+2,"solution error")

inferror = fhw01.inf_norm(e,1)
graph_f(inferror,nt,ot,dt,fnum+3,"solution error")


show()
