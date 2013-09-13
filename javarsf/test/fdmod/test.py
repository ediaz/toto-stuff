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

oz,dz,nz = 0.0,0.02,151
ox,dx,nx = 0.0,0.02,461
ot,nt,dt =    0.0, 1501, .002
fpeak = 13.0
ns = 81
###########################################################################
def main(args):
#  goMarmousi() # this creates the data and writes it to disk
  migshotsJava() # this one migrates the data


def migshotsJava():
  oz,dz,nz = 0.0,0.02,151
  ox,dx,nx = 0.0,0.02,461
  wav = zerofloat(nt,1)
  wav[0] = ricker()

  dat = readImage("allshots.dat",nt,nx,ns)
  src_coo = readImage("src_coo.dat",2,1,ns)
  rvr_coo = readImage("rvr_coo.dat",2,nx)

  vp2 = readImage('vinit.dat',nz,nx)

  mute = fdmod.Mute(ox,dx,nx,ot,dt,nt,3.2/fpeak,2*dz,1/fpeak)
  for isou in range(ns):
    xs = src_coo[isou][0][0]
    zs = src_coo[isou][0][1]
    mute.linear(xs,zs,1.6,dat[isou],dat[isou])

  rtm = fdmod.Rtm(ox,dx,nx,oz,dz,nz,ot,dt,nt,vp2,43,dat,ricker(),False)
  rtm.setSources(src_coo)   
  rtm.setReceivers(rvr_coo)
  rtm.setShotJump(20)
  rtm.goImage()
  img = rtm.getImage()
  plot(img,"")


def goMarmousi():
  marmfile='./vtrue.dat'
  mvel = readImage(marmfile,nz,nx)
  plotX(mvel,"marmousi")

  src_coo = zerofloat(2,1,ns)
  for isou in range(ns):
    xs = 1.*nx/ns*isou*dx +ox
    zs = 2*dz
    src_coo[isou][0][0] = xs 
    src_coo[isou][0][1] = zs
  wav = zerofloat(nt,1)
  wav[0] = ricker()

  rvr_coo = zerofloat(2,nx)
  for ix in range(nx):
    rvr_coo[ix][0] = ix*dx+ox
    rvr_coo[ix][1] = 2*dz
  fdobj = fdmod.Fdmod2d(ox,dx,nx,\
                        oz,dx,nz,\
                        ot,dt,nt,\
                        mvel,\
                        43,\
                        False) 

  os = 2*dx
  ds = nx/ns*dx
  fdobj.setReceivers(rvr_coo)

  dat = zerofloat(nt,nx,ns)
  sw = Stopwatch()
  sw.start()
  for isou in range(ns):
    sw.restart()
    xs = os+ds*isou
    zs = 2*dz
    fdobj.setSources(src_coo[isou],wav)
    fdobj.forwardAll()
    print "source %d out of %d done in %10.5g sec"%(isou,ns,sw.time())
    tmp =  fdobj.getData()
    xs = src_coo[isou][0][0]
    zs = src_coo[isou][0][1]
    copy(tmp,dat[isou])
  writeImage("allshots.dat",dat)
  writeImage("src_coo.dat",src_coo)
  writeImage("rvr_coo.dat",rvr_coo)

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
  fileName = name
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
  fileName = name
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
  pv.setPercentiles(2,99)
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
