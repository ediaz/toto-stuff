import helmhotzLibrary as lib
import numpy as np
import plot as plt
from copy import copy

'''
To run the figures in the report 
uncomment (remove #) in front of
the desired functions in function run()
'''
dirFig = 'report/Fig/'
def run():
  goImpulse(1,1,fname=dirFig+'FirstDerivative.pdf') # creates the impulse response for different Fornberg's filters for first derivatives
  goImpulse(2,2,fname=dirFig+'SecondDerivative.pdf') # creates the impulse response for different Fornberg's filters for second derivatives
#  goError(3)  # prints the tables for convergence analysis
#  goPluto(5)  # shows the Pluto related figures
  plt.show() 



def goPluto(fn=1):
  '''
  This function produces the response at 
  f = 10Hz for the Pluto synthetic model, which
  mimics the marine geology of the Gulf of Mexico.
  '''
  pgeo = {'o1':0,'d1':0.0608,'n1':125,
          'o2':0,'d2':0.0608,'n2':527}

  dwav = {'o1':0.0,'d1':0.004,'n1':251}

  
  pluto = lib.rsf2darray(**pgeo)  
  pluto.read('plutovel@')
  plt.plotvel(pluto,fn=fn,fname=dirFig+'Pluto.pdf')
  
  wav = lib.rsf1darray(**dwav)
  wav.ricker(15.0,0.15)
  plt.plot1d(wav,fn=fn+1,fname=dirFig+'twav.pdf')

  fftobj = lib.rsffft1(wav,Aunit='Hz',Aname='$f$')
  fir = fftobj.fft()
  fir.rsf = np.abs(fir.rsf)
  plt.plot1d(fir,fn=fn+2,fname=dirFig+'wwav.pdf')

  sx,sz = 16.,  0.05
  h_obj = lib.helmhotz(pluto,wav,sx,sz,nw=1,N=5,f=10.,free=True)

  u = h_obj.solve()
  
  wavefield = lib.rsf2darray(**pgeo)
  wavefield.put(u[0])
  plt.plotwavefield(wavefield,fn=fn+3,s=0.15,fname=dirFig+'free10Hz.pdf')

  h_obj = lib.helmhotz(pluto,wav,sx,sz,nw=1,N=5,f=10.,free=False)

  u = h_obj.solve()
  wavefield = lib.rsf2darray(**pgeo)
  wavefield.put(u[0])
  plt.plotwavefield(wavefield,fn=fn+10,s=0.15,fname=dirFig+'abc10Hz.pdf')



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
