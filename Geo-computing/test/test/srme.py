import sys
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from edu.mines.jtk.dsp import Conv
from edu.mines.jtk.dsp import FftReal 

############################################################################
#
#  look for Bill Dragoset's papers
#

def main(args):
  multiple()




def zoom2(o1f,n1f,o2f,n2f,x):

  y = zerofloat(n1f,n2f) 

  for i2 in range(n2f):
    for i1 in range(n1f):
      y[i2][i1] = x[i2+o2f][i1+o1f]

  return y




def multiple():
  n1,n2 = 2404,1001
  
  x = readImage()
  show2d(x,clip=0.1,title="x")

  y = autoConv(x)
  show2d(y,clip=0.1,title=" y=x*x")

  xx = zoom2(800,600,0,1001,x)
  yy = zoom2(800,600,0,1001,y)
  show2d(yy,clip=0.1,title="Zoom y=x*x")
  show2d(xx,clip=0.1,title="Zoom x")



def autoConv(x):

  n1,n2 = len(x[0]),len(x)
  y = zerofloat(n1,n2)

  for i2 in range(n2):
    Conv.conv(n1,0,x[i2],n1,0,x[i2],n1,0,y[i2]) 
    y[i2] = mul(-1.0,y[i2])
  return y








 
def like(x):
  n2 = len(x) 
  n1 = len(x[0])
  return zerofloat(n1,n2)

  






##########################################################################


def readImage():
  n1,n2 = 2404,1001
  x = zerofloat(n1,n2)
  ais = ArrayInputStream("noffset.dat",ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(x)
  ais.close()
  return x

def readMul():
  n1,n2 = 2404,1001
  x = zerofloat(n1,n2)
  ais = ArrayInputStream("noffset_ac.dat",ByteOrder.LITTLE_ENDIAN)
  ais.readFloats(x)
  ais.close()
  return x



def zoom2(o1f,n1f,o2f,n2f,x):

  y = zerofloat(n1f,n2f) 

  for i2 in range(n2f):
    for i1 in range(n1f):
      y[i2][i1] = x[i2+o2f][i1+o1f]

  return y

def plot(x,title):
  sp = SimplePlot.asPixels(x)
  sp.setTitle(title)
  sp.addColorBar()
  sp.setSize(800,1800)
  sp.plotPanel.setColorBarWidthMinimum(100)

def show2d(f,clip=1.0,title=None):
  print "show2d: f min =",min(f)," max =",max(f)


  clip *= (abs(min(f)) +abs(max(f)))*0.5


  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(f)
  if clip:
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(2,98)
  if title:
    sp.setTitle(title)
  sp.setSize(600,1100) 







#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
