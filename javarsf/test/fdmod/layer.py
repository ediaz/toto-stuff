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
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
import fdmod


'''
'''


plotWidth = 800
plotHeight = 800
plotPngDir = "./png/"
binDir = './layer/'


oz,dz,nz = 0.0,0.02,151
ox,dx,nx = 0.0,0.02,461
ot,nt,dt =    0.0, 1501, .002
fpeak = 11.0
ns = 1
###########################################################################
def main(args):
#  goMarmousi() # this creates the data and writes it to disk
#  migshotsJava() # this one migrates the data
  migshotsJava()


def goData():
  wav = zerofloat(nt,1)
  wav[0] = ricker()
  src_coo = zerofloat(2,1,ns)
  src_coo[0][0][0] = dx*nx/2.
  src_coo[0][0][1] = 2.*dz
  rvr_coo = zerofloat(2,nx)
  vel2 = getVel2()

  for ix in range(nx):
    rvr_coo[ix][0] = ix*dx+ox
    rvr_coo[ix][1] = 2*dz

  fdobj = fdmod.Fdmod2d(ox,dx,nx,oz,dx,nz,ot,dt,nt,vel2,43,False) 
  fdobj.setReceivers(rvr_coo)
  dat = zerofloat(nt,nx,ns)
  fdobj.setSources(src_coo[0],wav)
  fdobj.forwardAll()
  dat =  fdobj.getData()
  writeImage("allshots.dat",dat)
  writeImage("src_coo.dat",src_coo)
  writeImage("rvr_coo.dat",rvr_coo)

  plot(dat)
  

def getVel2():
  vel = zerofloat(nz,nx)
  for i2 in range(nx):
    for i1 in range(nz):
      if i1 < nz/2:
        vel[i2][i1] = 2.0
      else:
        vel[i2][i1] = 3.0
  return vel

def getVel():
  vel = zerofloat(nz,nx)
  for i2 in range(nx):
    for i1 in range(nz):
      vel[i2][i1] = 2.0
  return vel



def migshotsJava():
  oz,dz,nz = 0.0,0.02,151
  ox,dx,nx = 0.0,0.02,461
  wav = zerofloat(nt,1)
  wav[0] = ricker()

  dat = readImage("allshots.dat",nt,nx,ns)
  mask = readImage("mask.dat",nt,nx)
  src_coo = readImage("src_coo.dat",2,1,ns)
  rvr_coo = readImage("rvr_coo.dat",2,nx)

  vp2 = getVel2()

  mute = fdmod.Mute(ox,dx,nx,ot,dt,nt,10.4/fpeak,2*dz,1./fpeak)
  for isou in range(ns):
    xs = src_coo[isou][0][0]
    zs = src_coo[isou][0][1]
    mute.linear(xs,zs,3.,dat[isou],dat[isou])
    print len(dat),len(dat[0]),len(mask),len(mask[0])
    dat[isou]=mul(dat[isou],mask)
  plot(mask)
  plot(dat[isou]) 
  rtm = fdmod.Rtm(ox,dx,nx,oz,dz,nz,ot,dt,nt,vp2,43,dat,ricker(),False)
  rtm.setSources(src_coo)   
  rtm.setReceivers(rvr_coo)
  rtm.setShotJump(20)

  nh = 30
  jg = 1
  rtm.setExtendedImage(1, jg, (nx-nh)/jg, nh)
  rtm.goImage()
  img = rtm.getImage()
  ximg = rtm.getImageX()
  writeImage("ximage",ximg)

  plot(img,"")



def read_plotshots():
  shots = readImage("allshots.dat",nt,nx,ns)
  plot3d(readImage("allshots.dat",nt,nx,ns),0.1) 
  




########################
# Utility functions:   #
########################

def ricker():
  rick = zerofloat(nt)
  PI = acos(-1)
  delay = 1./fpeak
  s = (2*PI*PI*fpeak*fpeak)
  for i in range(nt):
    t = (ot +i*dt -delay)
    rick[i] = (1-s*t*t)*exp(-(PI*PI*fpeak*fpeak*t*t))
  return rick



def readImage(name,n1,n2,n3=None):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = binDir+name
  if n3==None:
    image = zerofloat(n1,n2)
  else:
    image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image


def writeImage(name,image):
  """ 
  Writes an image to a file with specified name.
  name: base name of image file; e.g., "tpgp"
  image: the image
  """
  fileName = binDir+name
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

def like(x):
  return zerofloat(len(x[0]),len(x))

###########################################################################
# Plotting

def plot2dmov(f,clip):
  min_clip= min(f)
  max_clip= max(f)
  for i3 in range(len(f)):
    plot(f[i3],"",None,min_clip,max_clip,clip)

def plotSequence(f,clip,png=None,title=""):
  n1 = len(f)
  s1 = Sampling(nt,dt,ot)
  p = panel()
  p.addTitle(title)
  p.setHLabel("t(s)")
  p.setHInterval(100*dt)
  p.setVLimits(-clip,clip)
  sv = p.addPoints(s1,f)
  frame(p,png)

def plot(x,title="",png=None,min_clip=None,max_clip=None,clip=1.):
  if min_clip==None:
    min_clip = min(x)
    max_clip = max(x)

  min_clip *= clip
  max_clip *= clip

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
  pv.setPercentiles(0,100)
#  pv.setClips(min_clip,max_clip)
  if not png==None:
    sp.paintToPng(720,3.33,plotPngDir+png+".png");



def plotX(x,title="",png=None,min_clip=None,max_clip=None):
  if min_clip==None:
    min_clip = min(x)
    max_clip = max(x)

  
  p = panel()
  sp = SimplePlot.asPixels(x)
  sp.setTitle(title)
  sp.setVLabel("z(km)")
  sp.setHLabel("x(km)")

  sp.addColorBar()
  sp.setSize(400*len(x)/len(x[0]),400)
  sp.plotPanel.setColorBarWidthMinimum(100)
  pv = sp.addPixels(x)
  pv.setColorModel(ColorMap.JET)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
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


def plot3d(f,clip=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])

  s1 = Sampling(n1,1,0)
  s2 = Sampling(n2,1,0)
  s3 = Sampling(n3,1,0)

  frame = SimpleFrame()
  ipg = frame.addImagePanels(s1,s2,s3,f)
  if clip:
    m  = max(abs(f))
    ipg.setClips(-clip*m,clip*m)
  ipg.setColorModel(ColorMap.GRAY)
  ipg.setSlices((n1/2),n2/2,n3/2)
  frame.orbitView.setScale(2.5)
  frame.setSize(1000,1000)
  if png:
    frame.paintToFile(plotPngDir+png+".png")

#########################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
