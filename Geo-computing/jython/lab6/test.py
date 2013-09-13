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

import lab6


'''
I took some plotting functions from Dave Hale's 
idh repository scripts (frame,panel,plot)
'''


plotWidth = 800
plotHeight = 800
plotPngDir = "./png/"
n1,dx,ox = 1001,1,-500
n2,dy,oy = 1001,1,-500

sigma=5

nrepeat = 8


###########################################################################
def main(args):
  #testGo1()
  #testGo2()
  #impulseResponses2D()

  goBenchmark() 
  


def testGo1():
  x = zerofloat(n1)
  x[n1/2] =1
  y = zerofloat(n1)
  y2 = zerofloat(n1)
  sigma=40
  gau = lab6.GaussianSmooth(sigma)
  gau.apply(x,y2)
  sequences = zerofloat(n1,7)
  i=0
  for nrepeat in [1,5,10,20,50,200]:
    exp = lab6.ExpSmooth(sigma,nrepeat)
    exp.apply(x,sequences[i])
    i+=1
  sequences[i]=y2
  plotSequences(sequences,png='Exp_Gauss')



def testGo2():
  x = randImage()
  y = like(x) 
  y2 = like(x) 

  exp = lab6.ExpSmooth(sigma,nrepeat)

  exp.apply(x,y)
  plot(y,"","Exp_2d_rand")

  gau = lab6.GaussianSmooth(sigma,False)
  gau.apply(x,y2)

  plot(y2,"","Gauss_2d_rand_0b")

  gau = lab6.GaussianSmooth(sigma,True)
  gau.apply(x,y2)
  plot(y2,"","Gauss_2d_rand_0sb")


def impulseResponses2D():
  y = zerofloat(n1,n2)
  finaly = like(y)
  sigma = 30
  nrepeat = 16

  centerX = [100,250,400]
  centerY = [100,250,400]
  repeat = [1,2,4,8,16,24,32,48,64]

  i=0
  # Impulse response of the 
  # Exponential filter for 
  # different number of smooth
  for i2 in centerY:
    for i1 in centerX:
      x=spike2d(i1,i2)
      exp = lab6.ExpSmooth(sigma,repeat[i])
      exp.apply(x,y)
  
      norm(y)
      add(y,finaly,finaly)
      i+=1

  clip = 1
  plot(finaly,"","Exp_2d_ir",0,clip)

  # Gaussian IR
  gau = lab6.GaussianSmooth(sigma,False)
  x = spike2d(250,250)
  gau.apply(x,y)
  plot(norm(y),"","Gauss_2d_ir_0b",0,clip)

 


 
def goBenchmark():
  global n1,n2
  n1,n2 = 1001,1001
  sigma = 10

  def benchmarkP(method,name,thread=1):
    x = randImage()
    sw = Stopwatch() #JTK stopwatch 
    sw.start()
    y = like(x)
    maxtime = 5.0
    nsmooth=0

    while (sw.time()<maxtime):
      method(x,y)    
      nsmooth += 1
    sw.stop()
    rate = nsmooth/sw.time()
    print ''+name+' benchmarking','rate= ',rate,'smooths/s'


  size=[100,200,500,1000]

  for n1 in size:
    n2=n1+1
    print "(n1,n2)=(%d,%d)"%(n1,n2)
    for nsmooth in [8,16,32]:
      e = lab6.ExpSmooth(sigma,nsmooth)
      benchmarkP(e.apply,'exponential n=%d'%nsmooth,1 )
    print '============================================================='
    for flen in [10,40,200]:
      e = lab6.GaussianSmooth(sigma,False)
      benchmarkP(e.apply,'Gaussian zb sigma=%d'%flen,1 )
      e = lab6.GaussianSmooth(sigma,True)
      benchmarkP(e.apply,'Gaussian zs sigma=%d'%flen,1 )

#<- End of goBenchmark function





################################
# Utility functions:

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


def plotSequence(f,dx=1,ox=0):
  n1 = len(f)

  s1 = Sampling(n1,dx,ox)
  p = panel()
  p.setHLabel("x")
  p.setHInterval(100*dx)
  p.setVLimits(min(f)*0.9,max(f)*1.1)
  sv = p.addPoints(s1,f)
  frame(p,None)




def plotLine(f,png=None):
  sv = panel.addPoints(1,0,f)
  sv.setLineColor(Color.BLACK)
  sv.setLineStyle(PointsView.Line.DOT)
  sv.setLineWidth(3)



def plot(x,title,png=None,min_clip=0.45,max_clip=0.55):
  p = panel() 
  sp = SimplePlot.asPixels(x)
  sp.setTitle(title)
#  sp.addColorBar()
  sp.setSize(800,800)
#  sp.plotPanel.setColorBarWidthMinimum(100)
  pv = sp.addPixels(x)
  pv.setClips(min_clip,max_clip)
  sp.paintToPng(720,3.33,plotPngDir+png+".png");


def randImage():
  x = zerofloat(n1,n2)
  rand = RandomFloat()
  for i2 in range(n2):
    for i1 in range(n1):
      x[i2][i1] = rand.uniform()
  return x


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
