from rsf.proj import *
from rsf.recipes import fdmod
import sys 

pt2in = 1./72
in2pt = 1./pt2in
four2three = {'width':356.96794*pt2in,
              'height':269.14662*pt2in
              }

sixteen2nine = {'width':448.01666*pt2in,
              'height':252.0748*pt2in
              }

parratio = {'4:3':four2three,
            '16:9':sixteen2nine}
            



class slideFig:
  def __init__(self,par,ratio='4:3',font=20,width=None,height=None,labelfat=2):
    self.par = par.copy()
    try:
      self.width = parratio[ratio]['width']
      self.height = parratio[ratio]['height']
      self.one20 = (parratio[ratio]['height']/font)*in2pt*0.5
    except:
      if width != None and height !=None:
        self.width = width
        self.height = height
        self.one20 = (height/font*in2pt*0.45)
      else:
        sys.exit("ERROR in slideFig class! \n either use default ratios 4:3 or 16:19, or you give me the dimensions")

    self.par['labelattr']='xll=2  parallel2=n labelsz=%f labelfat=%d titlesz=12 titlefat=3 '%(self.one20,labelfat)

  def printratio(self):
    print "width = %f in, height = %f in"%(self.width,self.height)

  def cgrey(self,custom,scalar=None):
    return self.hcgrey(custom,scalar)


  def hcgrey(self,custom,scalar=None):
    if scalar:
      self.par['height'] = self.par['ratio']*self.width*scalar
      if self.par['height'] >self.height : self.par['height']= self.height
    else:
      self.par['height'] = self.par['ratio']*self.width
      if self.par['height'] >self.height : self.par['height']= self.height
    return fdmod.cgrey(custom,self.par) 

  def vcgrey(self,custom,scalar=None):
    if scalar:
      self.par['height'] = self.par['ratio']*self.height*scalar
    else:
      self.par['height'] = self.par['ratio']*self.width
    if self.par['height'] >self.height : self.par['height']= self.height

    return fdmod.cgrey(custom,self.par) 




