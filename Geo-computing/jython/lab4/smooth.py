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

from lab4 import *

#############################################################################
def main(args):
  #goSmooth()
  #goSmooth1()
  #goSmooth2()
  #goDetrend()
  goBenchmark()

def goSmooth():
  x = randImage()
  a = 0.9
  y = smooth(a,x)
  plot(x,"input")
  plot(y,"smooth")

def goSmooth1():
  a = 0.9
  x = readImage()
  y1= smooth1(a,x)
  plot(x,"input")
  plot(y1,"smooth1")
  Dsp.smooth1P(a,x,y)
  plot(y,"smooth1P")

def goSmooth2():
  a = 0.9 
  x = readImage()
  y= smooth2S(a,x)
  plot(x,"input")
  plot(y,"smooth2")
  Dsp.smooth2P(a,x,y)
  plot(y,"smooth2P")

def goDetrend():
  a = 0.9
  x = readImage()
  y= detrend()
  plot(x,"input")
  plot(y,"detrend")


    
def goBenchmark():
  n1,n2 = 1001,1003
  x = randImage()
  a = 0.9

  def benchmarkP(method,name,thread=1):
      sw = Stopwatch() #JTK stopwatch 
      sw.start()
      nsmooth = 0
      y = like(x)
      maxtime = 3.0
      while (sw.time()<maxtime):
        method(a,x,y)    
        nsmooth += 1
      sw.stop()
      rate = (6.0e-6*n1*n2*nsmooth/sw.time())
      mean = Dsp.mean(y) 
      print ' '+name+': mflops=',int(rate),'mean=',mean
      return rate

  def benchmark2(methodS,nameS,methodP,nameP):
  
    print '========================================================='
    print 'benchmarking methods %s and %s'%(nameS,nameP)
    rs = benchmarkP(methodS,nameS)
    rp = benchmarkP(methodP,nameP,1 )
    print "rate= ",float(rp/rs)

  Dsp.thread=4
  benchmark2(Dsp.smooth1,'smooth1 ',Dsp.smooth1P,'smooth1P')
  benchmark2(Dsp.smooth2,'smooth2 ',Dsp.smooth2P,'smooth2P')
#<- End of goBenchmark function





#############################################################################

def like(x):
  return zerofloat(len(x[0]),len(x))

def detrend():
  a = 0.9
  x = readImage()
  y = smooth(a,x)
  return (sub(x,y))

def smooth(a,x):
  y = like(x)
  Dsp.smooth(a,x,y)
  return y

def smooth1(a,x):
  y = like(x)
  Dsp.smooth1(a,x,y)
  return y

def smooth1P(a,x):
  y = like(x)
  Dsp.smooth1P(a,x,y)
  return y

def smooth2(a,x):
  y = like(x)
  Dsp.smooth2(a,x,y)
  return y

def smooth2P(a,x):
  y = like(x)
  Dsp.smooth2P(a,x,y)
  return y

def readImage():
  n1,n2 = 750,600
  x = zerofloat(n1,n2)
  ais = ArrayInputStream("esteban.dat",ByteOrder.BIG_ENDIAN)
  ais.readFloats(x)
  ais.close()
  return x

def randImage():
  n1,n2 = 1001, 1003
  x = zerofloat(n1,n2)

  rand = RandomFloat()

  for i2 in range(n2):
    for i1 in range(n1):
      x[i2][i1] = rand.uniform()

  return x






def plot(x,title):
  sp = SimplePlot.asPixels(x)
  sp.setTitle(title)
  sp.addColorBar()
  sp.setSize(800,850)
  sp.plotPanel.setColorBarWidthMinimum(100)






#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
