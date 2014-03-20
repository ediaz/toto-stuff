import rsf.api as rsf
import numpy as np
from fdmod1dlib import *



par = rsf.Par()

Fvel = rsf.Input("vel")
Fden = rsf.Input("den")
Fwav = rsf.Input() # read wavelet from stdin
Fwave = rsf.Output() # output wavefield

ox = Fvel.float("o1")
dx = Fvel.float("d1")
nx = Fvel.int("n1")

ot = Fwav.float("o1")
dt = Fwav.float("d1")
nt = Fwav.int("n1")

source_x = par.float("sx",0.0) # source position
nb = par.int("nb",500) # boundary layer, is 1d u can go crazy

vel = np.zeros(nx,'f')
den = np.zeros(nx,'f') 

Fvel.read(vel)
Fden.read(den)

ns = 1 
source = np.zeros((ns,nt),'f')
Fwav.read(source[0,:])

sx = np.ones(ns)
sx[0] = source_x

wave = wave1d(ox,dx,nx,ot,dt,nt,vel,den,nb)
wave.init_source(source,sx)
movie = wave.apply() # output wavefield

Fwave.put("n2",nx)
Fwave.put("o2",ox)
Fwave.put("d2",dx)
Fwave.put("unit2",Fvel.string("unit1"))
Fwave.put("label2",Fvel.string("label1"))

Fwave.put("n1",nt)
Fwave.put("o1",ot)
Fwave.put("d1",dt)


Fwave.write(movie) 

