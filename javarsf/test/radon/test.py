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
from radon import *

'''
'''


plotWidth = 800
plotHeight = 800
plotPngDir = "./png/"

ox,dx,nx = -0.3,0.001,601
op,dp,np = -1.5,0.025,121
ot,nt,dt = 0.0, 1001, .004
fpeak = 13.0
ns = 81
###########################################################################
def main(args):
  data  = getslant()
  ssmod = zerofloat(nt,np)
  dataAdj= zerofloat(nt,nx)
  plot(data) 
  rad = Radon(ox,dx,nx,op,dp,np,ot,dt,nt) 
 
  ssmod = rad.Adjoint(data)
  plot(ssmod)
  dataAdj = rad.Forward(ssmod)
  plot(dataAdj)
  writeImage('adj.dat',ssmod)
  writeImage('dadj.dat',dataAdj)
def getslant():
  return readImage('slant.dat',nt,nx)

########################
# Utility functions:   #
########################




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
