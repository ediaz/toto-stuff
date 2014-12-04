import helmhotzLibrary as lib
import numpy as np
import plot as plt
from copy import copy
import matplotlib.pyplot as plt2
from mpl_toolkits.axes_grid1 import make_axes_locatable

'''
To run the figures in the report 
uncomment (remove #) in front of
the desired functions in function run()
'''
dirFig = 'report/Fig/'
def run():
  #goImpulse(1,1,fname=dirFig+'FirstDerivative.pdf') # creates the impulse response for different Fornberg's filters for first derivatives
  #goImpulse(2,2,fname=dirFig+'SecondDerivative.pdf') # creates the impulse response for different Fornberg's filters for second derivatives
#  goError(3)  # prints the tables for convergence analysis
  goPluto(5)  # shows the Pluto related figures
  #simple()
  plt.show() 


def simple():
  pgeo = {'o1':0,'d1':0.0608,'n1':15,
          'o2':0,'d2':0.0608,'n2':15}

  dwav = {'o1':0.0,'d1':0.004,'n1':251}

  vel  = np.zeros((15,15),'f') 
  vel  = 2.0
  
  vmod = lib.rsf2darray(**pgeo) 
  vmod.put(vel)

  wav = lib.rsf1darray(**dwav)
  wav.ricker(15.0,0.15)
  fftobj = lib.rsffft1(wav,Aunit='Hz',Aname='$f$')
  fir = fftobj.fft()
  fir.rsf = np.abs(fir.rsf)
  sx=0.12
  sz=0.02
  h_obj = lib.helmhotz(vmod,wav,sx,sz,nw=1,N=5,f=10.,free=True)
  h_obj.solve()
#
#  matrix = h_obj.lhs(13.0)
#
#  fig = plt.figure(fn,figsize=(10,5))
#  ax  = fig.add_subplot(111)
#  ax.imshow(matrix.todense().real)
#


def goPluto(fn=1):
  '''
  This function produces the response at 
  f = 10Hz for the Pluto synthetic model, which
  mimics the marine geology of the Gulf of Mexico.
  '''
  pgeo = {'o1':0,'d1':0.0608,'n1':15,
          'o2':0,'d2':0.0608,'n2':15}

  dwav = {'o1':0.0,'d1':0.004,'n1':251}

  
  pluto = lib.rsf2darray(**pgeo)  
  vel =  np.ones((pgeo['n2'],pgeo['n1']),'f')*2.0
  for i1 in range(pgeo['n2']):
    if i1>7:
      vel[i1,:] *= 1.5
  #pluto.read('plutovel@')
  pluto.put(vel)
  plt.plotvel(pluto,fn=fn)
  
  wav = lib.rsf1darray(**dwav)
  wav.ricker(15.0,0.15)
  plt.plot1d(wav,fn=fn+1)

  fftobj = lib.rsffft1(wav,Aunit='Hz',Aname='$f$')
  fir = fftobj.fft()
  fir.rsf = np.abs(fir.rsf)
  plt.plot1d(fir,fn=fn+2)

  par = pgeo.copy()
  sx,sz = par['n2']*0.5*par['d2'] , par['n1']*0.5*par['d1'] 
  h_obj = lib.helmhotz(pluto,wav,sx,sz,nw=1,N=11,f=10.,free=True)

  u = h_obj.solve()
  
  wavefield = lib.rsf2darray(**pgeo)
  wavefield.put(u[0])
  plt.plotwavefield(wavefield,fn=fn+3,s=0.15)

  h_obj = lib.helmhotz(pluto,wav,sx,sz,nw=1,N=5,f=10.,free=False)

  u = h_obj.solve()
  wavefield = lib.rsf2darray(**pgeo)
  wavefield.put(u[0])
  plt.plotwavefield(wavefield,fn=fn+10,s=0.15)

  A = h_obj.returnA()

  
  fig = plt2.figure(20,figsize=(10,5))
  ax  = fig.add_subplot(111)
  im = ax.imshow(A.todense().real)
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="2%", pad=0.05)
  cb = fig.colorbar(im,cax=cax)


def goConstant(order=4,fn=1):
  '''
  This function does the error analysis between 
  the computed solution and a known exact function.
  '''
  def dictconstant(h):
    mcons = {'o1':0.0,'d1':h,'n1':int(2./h)+1,
             'o2':0.0,'d2':h,'n2':int(2./h)+1}
    return mcons

  def uex(x,z,f=4.):
    return f*np.sin(f*np.pi*x)*np.sin(f*np.pi*z)

  def fex(x,z,w,f=4.):
    c = f*np.pi
    return (w*w-2*c*c)*np.sin(f*np.pi*x)*np.sin(f*np.pi*z)

  def error(diff):
    n = len(diff)
    in2 = 1./(n*n)
    return in2*np.sqrt((diff*diff).sum())

  def uofh(h):
    dic = dictconstant(h)
    oz,dz,nz = dic['o1'],dic['d1'],dic['n1']
    ox,dx,nx = dic['o2'],dic['d2'],dic['n2']
    x =  np.linspace(ox,1.,nx)
    z =  np.linspace(oz,1.,nz)
    x, z = np.meshgrid(x,z)
  
    f = 4.
    w = 0.
    ue = uex(x,z,f)
    vmin = ue.min(); vmax = ue.max()
    uexact = lib.rsf2darray(**dic)  
    uexact.put(ue)
    fexact = lib.rsf2darray(**dic)
    fexact.put(fex(x,z,w,f))
  
    vel = lib.rsf2darray(**dic)
    vel.put(np.ones((nx,nz)))
  
    hobj = lib.helmhotz(vel,None,0.5,0.5,nw=1,N=order+1,abc=False,verb=False)
  
    b = np.reshape(fexact.rsf,(nx*nz))
    u = hobj.solveE(w,b)
    sol = lib.rsf2darray(**dic)
    sol.put(u)
    return u -ue 
  
  def goPlot(h,fn):
    dic = dictconstant(h)
    oz,dz,nz = dic['o1'],dic['d1'],dic['n1']
    ox,dx,nx = dic['o2'],dic['d2'],dic['n2']
    x =  np.linspace(ox,1.,nx)
    z =  np.linspace(oz,1.,nz)
    x, z = np.meshgrid(x,z)
  
    f = 4.
    w = 0.
    ue = uex(x,z,f)
    vmin = ue.min(); vmax = ue.max()
    uexact = lib.rsf2darray(**dic)  
    uexact.put(ue)
    plt.plotvel(uexact,fn=fn,cbar=False,fname=dirFig+'known.pdf')

  Ks = [0,1,2,3]
  L2 = 1.

  for k in Ks:
    h = 0.04*order/pow(order,k)
    diff =  uofh(h)
    L2n = error(diff)
    print '%3d& %10.8f& %10.8f& %10.8g'%(int(2./h)+1,h,np.log(L2/L2n)/np.log(2),L2n)
    L2 = L2n
  goPlot(h,fn)
  

def goError(fn=1):
  i=0
  for o in [2,4]:
    print "convergence analysis for order=%d"%o
    print "  N&     h     &     r     &     e"
    goConstant(order=o,fn=fn+i)
    i+1  


def goImpulse(d=2,fn=1,morder=12,fname=None):
  diff = []*3
  order = range(2,morder+1,2)
  i=0  
  ir = lib.rsf1darray(0,1,200)
  for o in order:
    diff = lib.Diff(order=d,n=o+1,h=1.) 
    i+=1 
    a = np.zeros(201)
    a[100] = 1.
    ir.rsf = np.convolve(a,diff.filter(),mode='same')
    fftobj = lib.rsffft1(ir)
    fir = fftobj.fft()
    fir.rsf = np.abs(fir.rsf)
    plt.plot1d(fir,fn=fn)
  fir.rsf = (fir.grid()*2.*np.pi)**d
  plt.plot1d(fir,fn=fn,fname=fname)


#####################################################

run()
