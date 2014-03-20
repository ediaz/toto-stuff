import sys
import math
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

import lab7


plotWidth = 800
plotHeight = 800
plotPngDir = "./png/"

ox,nx,dx =    0.0, 401, 0.005
ot,nt,dt =    0.0, 1801, 0.001


###########################################################################
def main(args):
  for rho in [0.5]:
    awefd1d(rho) 

def awefd1d(rho):
  global plotPngDir

  tdelay = .05
  source,sourceX = d_dt_gauss_time(ox,tdelay,20)

  vel = zerofloat(nx)
  den = zerofloat(nx)

  for i in range(nx):
    vel[i] = 2.0
  for ix in range(nx/3):
    den[ix]=1
  for ix in range(nx/3,2*nx/3,1):
    den[ix]=2.
  for ix in range(2*nx/3,nx,1):
    den[ix]=3.

  fd1d = lab7.Awefd1dp(ox,dx,nx,ot,dt,nt,vel,den,source,sourceX,1,501)
  movie = fd1d.apply()
  plotSequence(source[0],max(abs(source[0])))
  plot(movie,"")

################################
# Utility functions:
def d_dt_gauss_time(_os,t0,f):
  sourceX = zerofloat(1)
  sourceX[0] = _os
  s = 1./(4*f)

  source = zerofloat(nt,1)
  is2 = 1./(2.*s*s)
  ct = 1.0/(math.sqrt(2*math.pi)*s)
  
  for it in range(nt):
    t = ot+it*dt
    source[0][it]  = math.exp(-(t-t0)*(t-t0)/(2*s*s))
    source[0][it] *= ct*(-2*(t-t0)/(2*s*s))
  return source,sourceX

def d_dt_gaussian(_os,ss,_ot,st):
  source = zerofloat(nx,nt)
  sourceX = zerofloat(nx)
  
  for isou in range(nx):
    x = isou*dx+ox
    sourceX[isou] = x

  cx = 1.0/(math.sqrt(2*math.pi)*ss)
  ct = 1.0/(math.sqrt(2*math.pi)*st)

  for it in range(nt):
    for isou in range(nx):
      x = isou*dx + ox
      t = ot+it*dt
      source[it][isou]  = math.exp(-(x-_os)*(x-_os)/(2*ss*ss))
      source[it][isou] *= ct*(-2*(t-_ot)/(2*st*st))
      source[it][isou] *= math.exp(-(t-_ot)*(t-_ot)/(2*st*st)) 

  return source,sourceX


def mean(x):
  n1= len(x[0])
  n2= len(x)
  m=0
  for i2 in range(n2):
    for i1 in range(n1):
      m+=x[i2][i1]
  return(m/(n1*n2))

def spike2d(i1,i2):
  x = zerofloat(n1,n2)
  x[i1][i2]=1
  return x

def norm(x):
  n2 = len(x)
  n1 = len(x[0])
  maxx=max(x)
  for i2 in range(n2):
    for i1 in range(n1):
      x[i2][i1]=x[i2][i1]/maxx
  return x

def like(x):
  return zerofloat(len(x[0]),len(x))

###########################################################################
# Plotting


def plotSequence(f,clip,png=None,title=""):
  n1 = len(f)
  s1 = Sampling(n1,dx,ox)
  p = panel()
  p.addTitle(title)
  p.setHLabel("x(km)")
  p.setHInterval(100*dx)
  p.setVLimits(-clip,clip)
  sv = p.addPoints(s1,f)
  frame(p,png)


def plot(x,title,png=None,min_clip=None,max_clip=None):
  if min_clip==None:
    min_clip = min(x)
    max_clip = max(x)

  p = panel()
  p.setHLabel("t(s)")
  p.setVLabel("x(km)")
 
  sp = SimplePlot.asPixels(x)
  sp.setTitle(title)
  sp.setVLabel("x(km)")
  sp.setHLabel("t(s)")

  sp.addColorBar()
  sp.setSize(800,800)
  sp.plotPanel.setColorBarWidthMinimum(100)
  pv = sp.addPixels(x)
  pv.setClips(min_clip,max_clip)
  if not png==None:
    sp.paintToPng(720,3.33,plotPngDir+png+".png");


def panel():
  p = PlotPanel(1,1,
    PlotPanel.Orientation.X1RIGHT_X2UP,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  return p

def frame(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSizeForSlide(1.5,1.2)
  frame.setSize(plotWidth,plotHeight)
  frame.setVisible(True)
  if png and plotPngDir:
    frame.paintToPng(720,3.3,plotPngDir+png+".png")
  return frame

def plotSequences(fs,png=None):
  nf = len(fs)
  n1 = len(fs[0])
  s1 = Sampling(n1,1.0,-(n1-1)/2.0)
  cs = ["r","g","b","c","m","k","y","b"]
  p  = panel()
  p.setHLabel("sample index")
  p.setHInterval(200.0)
  p.setVLimits(0.0,max(fs))
  for i in range(nf):
    sv = p.addPoints(s1,fs[i])
    if i >= len(cs):
      sv.setStyle(cs[(i-len(cs))%len(cs)]+'-')
    else:
      sv.setStyle(cs[i]+'-')
  frame(p,png)

#########################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
