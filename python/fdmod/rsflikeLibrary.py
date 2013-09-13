import numpy as np


class rsf2darray:
  def __init__(self,o1,d1,n1,o2,d2,n2,Aunit='km',\
                    Punit='km/s',prop='vel',name=''):
    self.name  = name # name of the file
    self.o1 = o1 # origin of the fast axis
    self.d1 = d1 # sampling interval
    self.n1 = n1 # number of samples of fast axis
    self.o2 = o2 # origin of the slow axis
    self.d2 = d2 # sampling interval of slow axis
    self.n2 = n2 # number of samples of slow axis
    self.Aunit = 'km' # axis unit
    self.Punit = 'km/s' # property unit (goes to the colorbar)
    self.prop = prop # property of array (default is vel)
    self.rsf = None # array to be read

  def read(self,fname):
    self.fname = fname # string variable with name to binary file
    self.rsf = np.fromfile(self.fname,'f')
    self.rsf = np.reshape(self.rsf,(self.n2,self.n1)) 

  def put(self,rsf):
    self.rsf = rsf   

  def geom(self):
    print 'o2=',self.o2,'d2=',self.d2,'n1=',self.n2
    print 'o1=',self.o1,'d1=',self.d1,'n1=',self.n1

class rsf1darray:
  def __init__(self,o1,d1,n1,Aunit='s',Aname='t'):
    self.o1 = o1 # origin of the fast axis
    self.d1 = d1 # sampling interval
    self.n1 = n1 # number of samples of fast axis
    self.rsf = None    
    self.m1 = o1+(n1-1)*d1
    self.Aunit = Aunit  
    self.Aname = Aname 
  def put(self,rsf):
    self.rsf = rsf
    if not len(rsf) == self.n1:
      print 'n1 attribute =! len(rsf)'

  def ricker(self,f,pt=0.0):
    t = self.grid()
    d = pt*np.ones(self.n1)
    s = 1./(4.*f)
    scale = 2./(np.sqrt(3*s)*np.pi**0.25)
    is2   = 1./(2*s*s)
    t2 = t-d
    st    = -t2*t2*2*is2
    expt  = np.exp(-t2*t2*is2)
    self.rsf   = (scale*(1+st)*expt)

  def geom(self):
    print 'o1=',self.o1,'d1=',self.d1,'n1=',self.n1

  def grid(self):
    return np.linspace(self.o1,self.m1,self.n1)


class rsffft1:
  '''
  input is a rsf1darray, it return an rsf1darray with the fft
  '''
  def __init__(self,rsf1d,Aunit='cycles/sample',Aname='$f_s$'):
    self.tarr = rsf1d
    self.ow = 0.0
    self.nyq = .5/rsf1d.d1
    self.Aunit = Aunit  
    self.Aname = Aname          
  def fft(self):
    from numpy.fft import rfft,irfft
    farr = rfft(self.tarr.rsf)
    nw = len(farr)
    dw = self.nyq/(nw-1)
    out = rsf1darray(self.ow,dw,nw,Aunit=self.Aunit,Aname=self.Aname)
    out.put(farr)
    return out  


