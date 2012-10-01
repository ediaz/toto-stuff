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

from lab3 import *

#############################################################################
def main(args):
  #goSmooth()
  #goSmooth1()
  #goSmooth2()
  #goSmooth2T()
  #goSmooth2S()
  #goDetrend()
  goBenchmark()

def goSmooth():
  a = 0.9
  x = readImage()
  y = smooth(a,x)
  plot(x,"input")
  plot(y,"smooth")

def goSmooth1():
  a = 0.9
  x = readImage()
  y1= smooth1(a,x)
  plot(x,"input")
  plot(y1,"smooth1")

def goSmooth2():
  a = 0.9 
  x = readImage()
  y1= smooth2S(a,x)
  plot(x,"input")
  plot(y1,"smooth2")


def goSmooth2T():
  a = 0.9 
  x = readImage()
  y1= smooth2T(a,x)
  plot(x,"input")
  plot(y1,"smooth2T")

def goSmooth2S():
  a = 0.9 
  x = readImage()
  y = smooth2S(a,x)
  plot(x,"input")
  plot(y,"smooth2S")


def goDetrend():
  a = 0.9
  x = readImage()
  y= detrend()
  plot(x,"input")
  plot(y,"detrend")
    
def goBenchmark():
  n1,n2 = 750,600
  x = readImage()
  a = 0.9

  def benchmarkP(method,name,thread=1):
      '''
      benchmark parallel code, if 
      thread=1 assumes it will benchmark
      sequential code
      '''
      Dsp.nthread =thread
      sw = Stopwatch() #JTK stopwatch 
      sw.start()
      nsmooth = 0
      y = like(x)
      maxtime = .001
      while (sw.time()<maxtime):
        method(a,x,y)    
        nsmooth +=1
      sw.stop()
      rate=int(6.0e-6*n1*n2*nsmooth/sw.time())
    
      print ' '+name+': mflops=','%4d'%rate,'mean=',Dsp.mean(y),'threads=',thread


  def benchmark2(methodS,nameS,methodP,nameP):
    '''
    Compare sequential implementation with parallel
    for different number of cores
    ''' 
    max_threads = 2 #Dsp.nthread
    print '========================================================='
    print 'benchmarking methods %s and %s'%(nameS,nameP)
    benchmarkP(methodS,nameS,1)
    for i in range(2,max_threads+1,1):
      benchmarkP(methodP,nameP,i )


  benchmark2(Dsp.smooth2,'smooth2 ',Dsp.smooth2P,'smooth2P')



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
  Dsp.smooth2(a,x,y)
  return y


def smooth2S(a,x):
  y = like(x)
  Dsp.smooth2S(a,x,y)
  return y

def smooth2SP(a,x):
  y = like(x)
  Dsp.smooth2SP(a,x,y)
  return y

def smooth2T(a,x):
  y = like(x)
  Dsp.smooth2T(a,x,y)
  return y

def smooth2TP(a,x):
  y = like(x)
  Dsp.smooth2TP(a,x,y)
  return y


def readImage():
  n1,n2 = 750,600
  x = zerofloat(n1,n2)
  ais = ArrayInputStream("esteban.dat",ByteOrder.BIG_ENDIAN)
  ais.readFloats(x)
  ais.close()
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
