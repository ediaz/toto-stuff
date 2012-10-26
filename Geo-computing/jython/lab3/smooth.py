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

##########################################################################
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

def goSmooth2T():
  a = 0.9 
  x = readImage()
  y= smooth2T(a,x)
  plot(x,"input")
  plot(y1,"smooth2T")
  Dsp.smooth2SP(a,x,y)
  plot(y,"smooth2TP")

def goSmooth2S():
  a = 0.9 
  x = readImage()
  y = smooth2S(a,x)
  plot(x,"input")
  plot(y,"smooth2S")
  Dsp.smooth2SP(a,x,y)
  plot(y,"smooth2SP")

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
      Dsp.nthread =thread

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

      print name+': mflops=','%4d'%int(rate),'mean=',mean,'threads=',thread
      return rate

  def benchmark2(methodS,nameS,methodP,nameP):
    max_threads = 4

    print '========================================================='
    print 'benchmarking methods %s and %s'%(nameS,nameP)
    
    for i in range(1,max_threads+1,1):
      rs = benchmarkP(methodS,nameS)
      rp = benchmarkP(methodP,nameP,i )
      print "rate= ",float(rp/rs)

  benchmark2(Dsp.smooth1,'smooth1 ',Dsp.smooth1P,'smooth1P')
  benchmark2(Dsp.smooth2T,'smooth2T ',Dsp.smooth2TP,'smooth2TP')

  #
  # parallelizing the transpose seems to improve the 
  # speed of the method
  #
  benchmark2(Dsp.smooth2T,'smooth2T ',Dsp.smooth2TPP,'smooth2TPP')
  benchmark2(Dsp.smooth2S,'smooth2S ',Dsp.smooth2SP,'smooth2SP')
  
  #
  # smooth2P seems to better than smooth2 for more than one thread,    
  # with 4 threads it slows down. Maybe, it just makes sense 
  # to use more threads if the image has a large n1.
  #
  # This is probably related with the overhead of creating and launching 
  # a thread with the run() method.
  #
  benchmark2(Dsp.smooth2,'smooth2 ',Dsp.smooth2P,'smooth2P')


#<- End of goBenchmark function





##########################################################################

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

########################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
