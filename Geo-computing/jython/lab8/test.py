import sys
import math

from java.awt import *
from java.io import *
from java.lang import *
from java.nio import *
from java.util import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
import lab8


plotWidth = 800
plotHeight = 800
plotPngDir = "./png/"
n1,d1,o1 = 251,1,-500
n2,d2,o2 = 252,1,-500
n3,d3,o3 = 253,1,-500

flag = -999.0

perc1d = 0.8
perc2d = 0.9995
perc3d = 0.9993
seed = 1678961
###########################################################################
def main(args):
#  interp1d()
#  interp2d()
#  interp3d()
#  dt3()
#  dt2()
#  dt1()
  galilee()


##########################################################################
#            Interpolation tests                                         #
##########################################################################
def interp3d():
  x = image2_3d()
  interp1 = zerofloat(n1,n2,n3)
  interp2 = zerofloat(n1,n2,n3)
  interp3 = zerofloat(n1,n2,n3)
  print "3D interpolation demo: perc=%f : %d out of n1*n2*n3=%d"%(perc3d*100,int((1-perc2d)*n1*n2*n3),n1*n2*n3)

  nn = lab8.NearestNeighbor(flag)
  nn.apply(x,interp1)
  plot3d(interp1,png="interp3d_nearest")
  print "Displaying NN"

  sb = lab8.SibsonInterpolation(flag)
#  sb.apply(x,interp2) 
#  print "Display Sibson"
#  plot3d(interp2)

  sb.applyP(x,interp3) 
  print "Display SibsonP"
  plot3d(interp3,png="interp3d_sibsonP")


def interp2d():
  x = test1_2d()
#  plot(x,"2D input")
  interp = like(x)

  print "2D interpolation demo: perc=%f : %d out of n1*n2=%d"%(perc2d*100,int((1-perc2d)*n1*n2),n1*n2)

  rx,x1,x2 = rsample2d(perc2d,flag,x)

  nn = lab8.NearestNeighbor(flag)
  nn.apply(rx,interp)
  plot(interp,"",png="interp2d_nearest",x1=x1,x2=x2)

  interp = like(x)
  sb = lab8.SibsonInterpolation(flag)
  sb.apply(rx,interp) 
  plot(interp,"",png="interp2d_sibson",x1=x1,x2=x2)

  interp2 = like(x)
  sb = lab8.SibsonInterpolation(flag)
  sb.applyP(rx,interp2) 
  plot(interp2,"",png="interp2d_sibsonP",x1=x1,x2=x2)





def interp1d():
  x = zerofloat(n1)
  y = zerofloat(n1)
  interp= zerofloat(n1)
  interp2= zerofloat(n1)
  r = RandomFloat(seed)
  s=n1/8.0
  for i1 in range(n1):
    y[i1] = exp(-(i1-n1/2)*(i1-n1/2)/(2*s*s))

  minval = min(y)
  maxval = max(y)

  plotSequence(x,minclip=minval,maxclip=maxval)
  for i in range(n1):
    if r.uniform() <perc1d:
      x[i] = flag
    else:
      x[i] = y[i]

  nn = lab8.NearestNeighbor(flag)
  nn.apply(x,interp)

  nn = lab8.SibsonInterpolation(flag)
  nn.apply(x,interp2)
  
  plotSequences2([y,x,x],png="sampled_function1d")

  plotSequences2([interp,interp2],png="interp_1d")

  i1 = zeroint(n1)
  d = zerofloat(n1)

  cpt = lab8.ClosestPointTransform(flag)
  cpt.apply(x,i1,d)
  plotSequences([d],png="distance_1d") 

##########################################################################
#                       Distance Transform Tests                         #
##########################################################################

def dt3():
  x = image2_3d()
  d = zerofloat(n1,n2,n3)
  d2 = zerofloat(n1,n2,n3)
  i1 = zeroint(n1,n2,n3)
  i2 = zeroint(n1,n2,n3)
  i3 = zeroint(n1,n2,n3)

  cpt = lab8.ClosestPointTransform(flag)
  sw =  Stopwatch()
  print "Closest point transform"
  print "n1=%d   n2=%d  n3=%d"%(n1,n2,n3)
  sw.start()
  cpt.apply(x,i1,i2,i3,d)
  print "time 3d Closest Point:",sw.time()
  plot3d(d)

  sw.reset()
  print "Closest point transform(naive)"
  print "n1=%d   n2=%d  n3=%d"%(n1,n2,n3)
  cpt_naive = lab8.Grid(flag)
  sw.start()
  cpt_naive.apply(x,i1,i2,i3,d2)
  print "time 3d Closest Point:",sw.time()

  plot3d(d2,png="distance3d_CPT")
  diff = sub(d2,d)

  print max(diff),min(diff)




def dt2():
  x = test2_2d()
  x,x1,x2 = rsample2d(perc2d,flag,x)
  d = like(x)
  d2 = like(x)
  i1 = zeroint(n1,n2)
  i2 = zeroint(n1,n2)


  cpt = lab8.ClosestPointTransform(flag)
  sw =  Stopwatch()
  print "Closest point transform"
  print "n1=%d   n2=%d "%(n1,n2)
  sw.start()
  cpt.apply(x,i1,i2,d)
  print "time 2d Closest Point:",sw.time()
  plot(d,"",png="distance2d_CPT")
  dc = like(x)

  for j2 in range(n2):
    for j1 in range(n1):
      dx = (i2[j2][j1]-j2)
      dz = (i1[j2][j1]-j1)
      dc[j2][j1] = sqrt(dx*dx*1.0 +dz*dz*1.0)

  cpt_naive = lab8.Grid(flag)
  cpt_naive.apply(x,i1,i2,d2)
  plot(d2,"",png="distance2d_naive")

  diff = like(x)

  diff = sub(d2,d)
  plot(diff,"DT naive - DT 2d")
  print "naive-dt2d",min(diff),max(diff)

  diff = sub(d2,dc)
  plot(diff,"DT naive - DT 2d cpt")
  print "naive - dt2d cpt ",min(diff),max(diff)

  print "2D interpolation demo: perc=%f : %d out of n1*n2=%d"%(perc2d*100,int((1-perc2d)*n1*n2),n1*n2)



def dt1():
  x = zerofloat(n1)
  d = zerofloat(n1)
  i1 = zeroint(n1)
  r = RandomFloat(seed)
 
  for i in range(n1):
    x[i] = i*(1.0 + 0.3*r.uniform())

  for i in range(n1):
    if r.uniform() <perc1d:
      x[i] = flag

  cpt = lab8.ClosestPointTransform(flag)
  sw =  Stopwatch()
  print "Closest point transform"
  print "n1=%d   "%(n1)
  sw.start()
  cpt.apply(x,i1,d)
  print "time 1d Closest Point:",sw.time()
  plotSequence(d)


def galilee():
  global flag

  flag =-9.0e9

  triplets = readImage()
  y = triplets[0]
  x = triplets[1]
  z = triplets[2]

  for i in range(len(z)):
    z[i]+= 212

  print "min x=",min(x),"max x=",max(x)
  print "min y=",min(y),"max y=",max(y)

  ox, oy = min(x), min(y)
  dx, dy = .02, .02 
  nx, ny = int((max(x)-ox)/dx)+1, int((max(y)-oy)/dy)+1
  print "ox=%g dx=%g nx=%d"%(ox,dx,nx)
  print "oy=%g dy=%g ny=%d"%(oy,dy,ny)

  image = grid(x,y,z,ox,dx,nx,oy,dy,ny)

  sx = Sampling(nx,dx,ox)
  sy = Sampling(ny,dy,oy)
  minc = min(z)
  maxc = max(z)

  xx=zerofloat(nx,ny)
  xx=add(-25,xx)

  plotGal(xx,"",sx,sy,minclip=minc,maxclip=maxc,nearest=True,\
          png="galilee_scatter",x1=x,x2=y)


  interp = zerofloat(nx,ny)
  nn = lab8.NearestNeighbor(flag)
  nn.apply(image,interp)
  plotGal(interp,"",sx,sy,minclip=minc,maxclip=maxc,nearest=True,\
          png="galilee_nearest")
  
  interp = zerofloat(nx,ny)
  si = lab8.SibsonInterpolation(flag)
  si.applyP(image,interp)
  plotGal(interp,"",sx,sy,minclip=minc,maxclip=maxc,png="galilee_sibson")
  plotGal(dfdx(interp),"d/dx Sibson Galilee bathymetry",sx,sy,-3,3)
  plotGal(dfdy(interp),"d/dy Sibson Galilee bathymetry",sx,sy,-3,3)





def dfdx(f):
  n2 = len(f)
  n1 = len(f[0])
  x = zerofloat(n1,n2)
  c2 = -1/12.0
  c1 = 2/3.0 
  for i2 in range(2,n2-2,1):
    for i1 in range(n1):
      x[i2][i1] = c2*(f[i2+2][i1]-f[i2-2][i1])+c1*(f[i2+1][i1]-f[i2-1][i1])

  return x

def dfdy(f):
  n2 = len(f)
  n1 = len(f[0])
  x = zerofloat(n1,n2)
  c2 = -1/12.0
  c1 = 2/3.0 
  for i2 in range(n2):
    for i1 in range(2,n1-2,1):
      x[i2][i1] = c2*( f[i2][i1+2]-f[i2][i1-2])+c1*(f[i2][i1+1]-f[i2][i1-1])

  return x



def grid(x,y,z,ox,dx,nx,oy,dy,ny):
  nk = len(x)
  n1 = nx
  n2 = ny
  nn = zeroint(n1,n2)
  image = zerofloat(n1,n2)

  for n in range(nk):
    xi = int(((x[n]-ox)/dx)) 
    yi = int(((y[n]-oy)/dy)) 
      
    image[yi][xi] += z[n]  
    nn[yi][xi] += 1

  for i2 in range(n2):
    for i1 in range(n1):
      if nn[i2][i1] != 0:
        image[i2][i1] = image[i2][i1]/nn[i2][i1]
      else:
        image[i2][i1] = flag
  return image



################################
#       test images            #
################################
def image1_3d():
  x = zerofloat(n1,n2,n3)
  r = RandomFloat(seed)

  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        c3 = i3-n3/2.
        c2 = i2-n2/2.
        c1 = i1-n1/2.
        d = 1.0*sqrt(c3*c3+c2*c2+c1*c1)
        if d<n1/4 :
          x[i3][i2][i1] = flag
        else:
          x[i3][i2][i1] = d*(1 +.5*r.uniform())
  return x




def image2_3d():
  x = zerofloat(n1,n2,n3)
  r = RandomFloat(seed)
  s1 = n1/4
  s2 = n2/4
  s3 = n3/4
  c3 = n3/2.
  c2 = n2/2.
  c1 = n1/2.

  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        if r.uniform() < perc3d :
          x[i3][i2][i1] = flag
        else:
          x3 = (i3-c3)
          x2 = (i2-c2)
          x1 = (i1-c1)
          e= (x3*x3)/(2*s3*s3)+(x2*x2)/(2*s2*s2)+(x1*x1)/(2*s1*s1)
          x[i3][i2][i1] = exp(-e)

  return x




def test1_3d():
  x = zerofloat(n1,n2,n3)
  i3cycles = 1.0/n3*2.0*PI
  i2cycles = 1.0/n2*2.0*PI
  i1cycles = 1.0/n1*2.0*PI

  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        x[i3][i2][i1] =  sin(i3cycles*i3)*sin(i2cycles*i2)*sin(i1cycles*i1)  
  return x



def test1_2d():
  x = zerofloat(n1,n2)
  i2cycles = 1.0
  i1cycles = 1.0
  for i2 in range(n2):
    for i1 in range(n1):
      x[i2][i1] = sin((i2cycles/n2)*2.0*PI*i2)*sin((i1cycles/n1)*2.0*PI*i1)
  return x

def test2_2d():
  x = zerofloat(n1,n2)
  s1 = n1/4
  s2 = n2/4
  for i2 in range(n2):
    for i1 in range(n1):
      x[i2][i1] = exp(-((i2-n2/2.0)*(i2-n2/2.)/(2.*s2*s2)+(i1-n1/2.)*(i1-n1/2.)/(2.*s1*s1)))
  return x


################################
# Utility functions:
################################


def rsample2d(perc,flag,x):
  y = like(x)
  y = x
  #perc : percentage of sample to discard
  if perc>1: perc=1
  nr = int(n1*n2*perc)
  sample = range(n1*n2)

  rand = RandomFloat()
  n=n1*n2
  ik = 0
  for i2 in range(n2):
    for i1 in range(n1):
      if(rand.uniform()<perc):
        y[i2][i1] = flag
      else:
        ik+=1
  
  x1 = zerofloat(ik)
  x2 = zerofloat(ik)
  ik=0
  for i2 in range(n2):
    for i1 in range(n1):
      if(y[i2][i1] != flag):
        x1[ik] = i1
        x2[ik] = i2
        ik+=1


  return y,x1,x2

def like(x):
  return zerofloat(len(x[0]),len(x))



def readImage():
  n1t,n2t = 130646,3
  x = zerofloat(n1t,n2t)
  ais = ArrayInputStream("galilee.dat",ByteOrder.BIG_ENDIAN)
  ais.readFloats(x)
  ais.close()
  return x

######################################################################
# Plotting

def plotSequence(f,dx=1,ox=0,minclip=None,maxclip=None):
  n1 = len(f)

  s1 = Sampling(n1,dx,ox)
  p = panel()
  p.setHLabel("x")
  p.setHInterval(100*dx)
  if minclip==None:
    minclip = min(f)
    maxclip = max(f)

  p.setVLimits(minclip,maxclip)
  sv = p.addPoints(s1,f)
  frame(p,None)




def plotLine(f,png=None):
  sv = panel.addPoints(1,0,f)
  sv.setLineColor(Color.BLACK)
  sv.setLineStyle(PointsView.Line.DOT)
  sv.setLineWidth(3)



def plot(x,title,minclip=None,maxclip=None,gray=False,nearest=True,png=None,x1=None,x2=None):
  if minclip==None:
    minclip=min(x)
    maxclip=max(x)
    if(min(x) == max(x)):
      minclip=-0.0001
      maxclip=+0.0001
  

  p = panel() 
  sp = SimplePlot.asPixels(x)
  sp.setTitle(title)
  sp.addColorBar()
  sp.setSize(800,800)
  sp.plotPanel.setColorBarWidthMinimum(100)
  pv = sp.addPixels(x)
  if gray:
    pv.setColorModel(ColorMap.GRAY)
  else:
    pv.setColorModel(ColorMap.JET)
  pv.setClips(minclip,maxclip)
  if nearest:
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if x1 and x2:
    pntv = sp.addPoints(x1,x2)
    pntv.setLineStyle(PointsView.Line.NONE);
    pntv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pntv.setMarkSize(5.0);

  if png:
    sp.paintToPng(720,3.33,plotPngDir+png+".png");



def plotGal(x,title,sx,sy,minclip=None,maxclip=None,nearest=False,png=None,x1=None,x2=None):
  if minclip==None:
    minclip=min(x)
    maxclip=max(x)

  p = panel() 
  sp = SimplePlot.asPixels(sx,sy,x)
  sp.setTitle(title)
  sp.addColorBar()
  sp.setSize(800,800)
  sp.plotPanel.setColorBarWidthMinimum(100)
  pv = sp.addPixels(sx,sy,x)
#  pv.setColorModel(ColorMap.GRAY)
  pv.setColorModel(ColorMap.JET)
  pv.setClips(minclip,maxclip)
#  sp.setVRotated(True)
  nearest=False
  if nearest:
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  sp.setHLabel("East (km)")
  sp.setVLabel("North (km)")
  if x1 and x2:
    pntv = sp.addPoints(x1,x2)
    pntv.setLineStyle(PointsView.Line.NONE);
    pntv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pntv.setMarkSize(1.);


  if png:
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


def plotSequences2(fs,png=None):
  nf = len(fs)
  n1 = len(fs[0])
  s1 = Sampling(n1,1.0,-(n1-1)/2.0)
  cs = [Color.BLUE,Color.RED,Color.GREEN,Color.BLACK]
  p  = panel()
  p.setHLabel("sample index")
  p.setHInterval(200.0)
  p.setVLimits(0.0,max(fs))
  print nf
  for i in range(nf):
    sv = p.addSequence(s1,fs[i])
    sv.setColor(cs[i])
  frame(p,png)



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
    ipg.setClips(-clip,clip)
  ipg.setColorModel(ColorMap.JET)
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
