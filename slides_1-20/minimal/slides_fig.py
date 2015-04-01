from rsf.proj import *
from rsf.recipes import fdmod
import sys 

pt2in = 1./72
in2pt = 1./pt2in
four2three = {'width':356.96794*pt2in,
              'height':269.14662*pt2in
              }

sixteen2nine = {'width':455.24408*pt2in,
              'height':256.0748*pt2in
              }

parratio = {'4:3':four2three,
            '16:9':sixteen2nine}
            



class slideFig:
  def __init__(self,par,ratio='4:3',font=20,width=None,height=None,labelfat=2):
    self.par = par.copy()
    try:
      self.width = parratio[ratio]['width']
      self.height = parratio[ratio]['height']
      self.one20 = (parratio[ratio]['height']/font)*in2pt*0.5*0.7
    except:
      if width != None and height !=None:
        self.width = width
        self.height = height
        self.one20 = (height/font*in2pt*0.45)
      else:
        sys.exit("ERROR in slideFig class! \n either use"+ 
                 " default ratios 4:3 or 16:19, or you give me the dimensions")

    self.par['labelattr']='xll=2  parallel2=n labelsz=%f labelfat=%d titlesz=12 titlefat=3 '%(self.one20,labelfat)

  def print_dimensions(self):
    '''
    This function prints the dimensions of the selected slide
    '''
    print "width = %f in, height = %f in"%(self.width,self.height)



  ###################################################
  #
  #      Interface to fdmod.cgrey()
  #
  ###################################################

  def cgrey(self,custom,scalar=None,horizontal=True):
    '''
    This is equivalent to fdmod.cgrey('',par) 
    scalar:  horizontal scalar for the figure space
    '''
    self._scale(scalar,horizontal)
    return fdmod.cgrey(custom,self.par)

  def hcgrey(self,custom,scalar=None):
    '''
    read cgrey
    h stands for horizontal scalar
    '''
    self._h_scale(scalar)
    return self.cgrey(custom,scalar,True)

  def vcgrey(self,custom,scalar=None):
    '''
    read cgrey
    v stands for vertical scalar
    '''
    self._v_scale(scalar)
    return self.cgrey(custom,scalar,False)



  ###################################################
  #
  #      Interface to fdmod.dgrey()
  #
  ###################################################

  def dgrey(self,custom,scalar=None,horizontal=True):
    '''
    This is equivalent to fdmod.cgrey('',par) 
    scalar:  horizontal scalar for the figure space
    '''
    self._scale(scalar,horizontal)
    toplot = ' plotfat=1 screenratio=%(ratio)g screenht=%(height)g '%self.par +custom
    return fdmod.dgrey(toplot,self.par)

  def hdgrey(self,custom,scalar=None):
    '''
    read cgrey
    h stands for horizontal scalar
    '''
    self._h_scale(scalar)
    return self.dgrey(custom,scalar,True)

  def vdgrey(self,custom,scalar=None):
    '''
    read cgrey
    v stands for vertical scalar
    '''
    self._v_scale(scalar)
    return self.dgrey(custom,scalar,False)


  ###################################################
  #
  #      Interface to fdmod.egrey()
  #
  ###################################################

  def egrey(self,custom,scalar=None,horizontal=True):
    '''
    This is equivalent to fdmod.cgrey('',par) 
    scalar:  horizontal scalar for the figure space
    '''
    self._scale(scalar,horizontal)
    toplot = ' plotfat=1 screenratio=%(ratio)g screenht=%(height)g '%self.par +custom
    return fdmod.egrey(toplot,self.par)

  def hegrey(self,custom,scalar=None):
    '''
    read cgrey
    h stands for horizontal scalar
    '''
    self._h_scale(scalar)
    return self.egrey(custom,scalar,True)

  def vegrey(self,custom,scalar=None):
    '''
    read cgrey
    v stands for vertical scalar
    '''
    self._v_scale(scalar)
    return self.egrey(custom,scalar,False)



  ###################################################
  #
  #      Interface to fdmod.waveplot()
  #
  ###################################################

  def waveplot(self,custom,scalar=None,horizontal=True):
    '''
    This is equivalent to fdmod.cgrey('',par) 
    scalar:  horizontal scalar for the figure space
    '''
    self._scale(scalar,horizontal)
    toplot = ' plotfat=1 screenratio=%(ratio)g screenht=%(height)g '%self.par +custom
    return fdmod.waveplot(toplot,self.par)

  def hwaveplot(self,custom,scalar=None):
    '''
    read cgrey
    h stands for horizontal scalar
    '''
    self._h_scale(scalar)
    return self.waveplot(custom,scalar,True)

  def vwaveplot(self,custom,scalar=None):
    '''
    read cgrey
    v stands for vertical scalar
    '''
    self._v_scale(scalar)
    return self.waveplot(custom,scalar,False)

  ###################################################
  #
  #      Interface to fdmod.spectrum()
  #
  ###################################################

  def spectrum(self,custom,scalar=None,horizontal=True):
    '''
    This is equivalent to fdmod.cgrey('',par) 
    scalar:  horizontal scalar for the figure space
    '''
    self._scale(scalar,horizontal)
    return fdmod.spectrum(custom,self.par)

  def hspectrum(self,custom,scalar=None):
    '''
    read cgrey
    h stands for horizontal scalar
    '''
    self._h_scale(scalar)
    return self.spectrum(custom,scalar,True)

  def vspectrum(self,custom,scalar=None):
    '''
    read cgrey
    v stands for vertical scalar
    '''
    self._v_scale(scalar)
    return self.spectrum(custom,scalar,False)



  ###################################################
  #
  #      Interface to fdmod.ssplot and fdmod.rrplot
  #
  ###################################################


  def ssplot(self,custom,scalar,horizontal=True):
    self._scale(scalar,horizontal)
    return fdmod.ssplot(custom,self.par)
    
  def hssplot(self,custom,scalar):
    return self.ssplot(custom,scalar,True)

  def vssplot(self,custom,scalar):
    return self.ssplot(custom,scalar,False)
    
  def rrplot(self,custom,scalar,horizontal=True):
    self.scale(horizontal)
    return fdmod.rrplot(custom,self.par)
    
  def hrrplot(self,custom,scalar):
    return self.rrplot(custom,scalar,True)

  def vrrplot(self,custom,scalar):
    return self.rrplot(custom,scalar,False)
    


  def get_par(self):
    '''
    returns the current state of par file
    '''
    return self.par

  ## Private functions user should not call these:
  def _scale(self,scalar,horizontal=True):
    if horizontal:
      self._h_scale(scalar)
    else:
      self._v_scale(scalar)

    

  def _h_scale(self,scalar):
    if scalar:
      self.par['height'] = self.par['ratio']*self.width*scalar
      if self.par['height'] >self.height : self.par['height']= self.height
    else:
      self.par['height'] = self.par['ratio']*self.width
      if self.par['height'] >self.height : self.par['height']= self.height
    
  def _v_scale(self,scalar):
    if scalar:
      self.par['height'] = self.height*scalar
    else:
      self.par['height'] = self.par['ratio']*self.width
    if self.par['height'] >self.height : self.par['height']= self.height
    

