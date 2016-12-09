try: from rsf.cluster import *
except: from rsf.proj import *
from rsf.recipes import fdmod,mpi
from aweOP import aweop2d
from eicOP import eicop2d 
from cicOP import cicop2d
import ciputils
from injOP import injop2d




class wetID:
  def __init__(self,dS,ss,dR,rr,gg,par,mask=None,customfd='nb=20 order=15 jsnap=1'):
    self.par = par.copy()
 
    self.dS = dS
    self.ss = ss

    self.dR = dR
    self.rr = rr
  
    self.nexp = len(ss)
    self.mask = mask

    self.gg = gg 

    self.setSaveDir()
    self.setPrec()
    self.cip = False 
    self.setImPower()
    self.setProb()
    self.customfd = customfd

  def setPrec(self,rect1=1,rect2=1,repeat=1):
    self.prec = 'smooth rect1=%d rect2=%d repeat=%d'%(rect1,rect2,repeat)

  def setPicker(self,parcip):
    self.cip = True
    self.parcip = parcip.copy()

  def setCompute(self,nodes=5,time=1):
    self.nodes = nodes
    self.ipn =  int(self.nexp/nodes)
    self.time = time
    self.cluster = True

  def setImPower(self,im=.1):
    self.im = im



  def setSaveDir(self,sdir='images/'):
    self.sdir = sdir
  
  def stateVariables(self,m):
    self.tag(m)

    images = [] 
    eimages = [] 
    Fork(time=self.time,ipn=self.ipn,nodes=self.nodes)
    for iexp in range(self.nexp):
      stag = '-%03d'%iexp
      ss   = self.ss[iexp] 
      rr   = self.rr[iexp]
      wav  = self.dS[iexp]
      dat  = self.dR[iexp]
      vel  = 'tmp/vel-'+self.mtag+stag 
      swfl = 'tmp/swfl-'+self.mtag+stag
      rwfl = 'tmp/rwfl-'+self.mtag+stag
      img  = 'tmp/img-'+self.mtag+stag

      self.getVel(iexp,m,vel)

      Lswfl = aweop2d(vel,self.par,ss,'','',custom=self.customfd)
      Lswfl.FORW(wav,swfl)
    
      Lrwfl = aweop2d(vel,self.par,'',rr,'',custom=self.customfd)
      Lrwfl.ADJT(rwfl,dat)
    
      Rop = cicop2d(swfl,self.par)
      Rop.FORW(rwfl,img)
      images.append(img)
      Iterate()
    Join()


    image = self.sdir+'Image-'+self.mtag
    mpi.gridandstack(image,images,
     8,'tmp/img-'+self.mtag+'-%03d',
     self.par['nz'],self.par['oz'],self.par['dz'],
     self.par['nx'],self.par['ox'],self.par['dx'],
     1.0         ,0.0         , 0.1,
     nf=self.par['ns'],of=0,jf=1)

    #Flow(image,images,'add ${SOURCES[1:-1]}')
    self.image= image

    if self.cip:
      gg = self.sdir+'/gg-'+self.mtag
      self.getLocations(gg,image)
      self.gg = gg
    else:
      gg = self.gg

    Fork(time=self.time,ipn=self.ipn,nodes=self.nodes )
    for iexp in range(self.nexp):
      stag = '-%03d'%iexp
      swfl = 'tmp/swfl-'+self.mtag+stag
      rwfl = 'tmp/rwfl-'+self.mtag+stag
      eimg = 'tmp/eic-'+self.mtag+stag

      ERop = eicop2d(swfl,gg,self.par)
      ERop.FORW(rwfl,eimg)
      eimages.append(eimg)
      Iterate()
    Join()
    self.eimage = self.sdir+'EImage-'+self.mtag
    Flow(self.eimage,eimages,'add ${SOURCES[1:-1]}')

  def setMute(self,mute=None):
    self.mute = mute

  def computePenalty(self,m):
    self.penalty='tmp/penalty'+self.mtag
      
    Flow(self.penalty,[m,self.gg],
      '''
      python ./getpenalty.py nht=%(nht)d nhx=%(nhx)d nhz=%(nhz)d 
                             dht=%(dt)g  dhx=%(dx)g  dhz=%(dz)g
             cip=${SOURCES[1]}
      '''%self.par)

    self.residual = 'images/res-'+self.mtag
    self.res = 'images/resx-'+self.mtag
    sources = [self.eimage,self.penalty]

    Flow([self.residual,self.res],sources,
      'add mode=p ${SOURCES[1]}> ${TARGETS[1]} && add mode=p ${SOURCES[1]} <${TARGETS[1]}')

#    convert = '''
#    window |put n1=%d d1=%g o1=%g n2=%d d2=%g o2=%g n3=%d d3=%g o3=%g|transp
#    '''%(self.par['nhx']*2+1,self.par['dx'],-self.par['nhx']*self.par['dx'],
#         self.par['nz'],self.par['dz'],self.par['oz'],
#         self.par['ncig'],self.par['dcig'],self.par['ocig'])
#
#    convertBack = '''
#    transp |window |put n3=1 n2=1 n4=%d |transp 
#    '''%(self.par['nz']*self.par['ncig'])
#
#    Flow([self.residual,self.res],sources,
#      convert+
#      '''|math output="abs(x2)*input"|slant dp=0.01 p0=-3.0 np=321 adj=y >${TARGETS[1]} &&
#       sfslant adj=n x0=%g dx=%g nx=%d <${TARGETS[1]} |math output="abs(x2)*input"|'''%
#      (self.par['dx']*self.par['nhx']*-2.0,self.par['dx']*2,self.par['nhx']*2+1)+
#      convertBack)
#
  def setIllu(self,illu=False):
    self.doIllu = illu

  def computeIlluPenalty(self,m):
    self.penalty=self.sdir+'penalty'+self.mtag
    gg = self.gg
    
    PSF = 'tmp/PSFImage-'+self.mtag
    Flow(PSF,[m,self.gg],
      '''
      python ./getpenalty.py nht=%(nht)d nhx=%(nhx)d nhz=%(nhz)d 
                             dht=%(dt)g  dhx=%(dx)g  dhz=%(dz)g
             cip=${SOURCES[1]} spike=y
      '''%self.par) 

    Fork(time=self.time,ipn=self.ipn,nodes=self.nodes )
    for iexp in range(self.nexp):
      stag = '-%03d'%iexp
      rr   = self.rr[iexp]
      swfl = 'tmp/swfl-'+self.mtag+stag
      arillu = 'tmp/arillu-'+self.mtag + stag
      arilluwfl = 'tmp/arilluwfl-'+self.mtag + stag
      mdata = 'tmp/data-rillu-'+self.mtag + stag
      fmdata = 'tmp/fdata-rillu-'+self.mtag + stag
      vel  = 'tmp/vel-'+self.mtag+stag 
      rillu= 'tmp/rillu-'+self.mtag+stag
      gg   = self.gg

      ERop = eicop2d(swfl,gg,self.par)
      ERop.ADJT(arillu,PSF)

      Lswfl = aweop2d(vel,self.par,'',rr,'',custom=self.customfd)
      Lswfl.FORW(arillu,mdata)
      Lswfl.ADJT(arilluwfl,mdata)

      ERop.FORW(arilluwfl,rillu)
      Iterate()
    Join()

    Flow('images/rilluall-'+self.mtag,
         map(lambda x:'tmp/rillu-'+self.mtag+'-%03d'%x,range(self.nexp)),
         '''
         add ${SOURCES[1:%d]} 
         '''%(self.nexp))

    if self.minimize:
      Flow(self.penalty,'images/rilluall-'+self.mtag,
          '''
           transp plane=12|
         envelope |
         transp plane=12|
         transp plane=31|
         envelope |
         transp plane=31|
         smooth rect1=4 rect2=4 rect3=4|
         scale axis=3|
         math output="(1-input)"|
         math output="input*sqrt(x2^2+(3*x3)^2)" 
         ''')
    else:
      Flow(self.penalty,'images/rilluall-'+self.mtag,
          '''
           transp plane=12|
         envelope |
         transp plane=12|
         transp plane=31|
         envelope |
         transp plane=31|
         smooth rect1=5 rect2=5 rect3=5|
         scale axis=3
         ''')


    self.residual = 'images/res-'+self.mtag
    self.res = 'images/resx-'+self.mtag
    sources = [self.eimage,self.penalty]

    Flow([self.residual,self.res],sources,
      'add mode=p ${SOURCES[1]}> ${TARGETS[1]} && add mode=p ${SOURCES[1]} <${TARGETS[1]}')

  def setProb(self,min=True):
    self.minimize=min


  def of(self,OF,m):
    self.tag(m)

    self.stateVariables(m)
    if self.doIllu:
      self.computeIlluPenalty(m)
    else:
      self.computePenalty(m)

    
    #Flow(OF,[self.penalty,self.eimage],
    #  'add mode=p ${SOURCES[1]} |math output="input^2"|stack axis=0 norm=n|scale rscale=%g'%self.scalar)
    addition = ''
    if not self.minimize: '|scale rscale=-1' 
    Flow(OF,[self.res],
       'add mode=p $SOURCE |stack axis=0 norm=n |scale rscale=%g'%self.scalar+addition)


  def adjointStateVariables(self,m,stateTag=None):
    if stateTag==None: stateTag = self.mtag
    sgrads  = [] 
    rgrads  = [] 
    sgradsIm  = [] 
    rgradsIm  = [] 
    Fork(time=self.time,ipn=self.ipn,nodes=self.nodes )
    for iexp in range(self.nexp):
      stag = '-%03d'%iexp
      ss   = self.ss[iexp] 
      rr   = self.rr[iexp]
      gg   = self.gg

      vel  = 'tmp/vel-'+self.mtag+stag 

      res  = self.residual 

      # forward variables
      swfl = 'tmp/swfl-'+stateTag+stag
      rwfl = 'tmp/rwfl-'+stateTag+stag

      # adjoint sources 
      adjS = 'tmp/adj-swfl-'+self.mtag+stag
      adjSim = 'tmp/adj-imswfl-'+self.mtag+stag
      adjR = 'tmp/adj-rwfl-'+self.mtag+stag
      adjRim = 'tmp/adj-imrwfl-'+self.mtag+stag


      # adjoint state variables
      aswfl = 'tmp/aswfl-'+self.mtag+stag
      asimp = 'tmp/asimwfl-'+self.mtag+stag
      arwfl = 'tmp/arwfl-'+self.mtag+stag
      arimp = 'tmp/arimwfl-'+self.mtag+stag

      # gradients
      sgrad = 'tmp/sgrad-'+self.mtag+stag
      rgrad = 'tmp/rgrad-'+self.mtag+stag
      sgradIm = 'tmp/sgradIm-'+self.mtag+stag
      rgradIm = 'tmp/rgradIm-'+self.mtag+stag
    
      Timage = self.sdir+'Image-'+self.mtag
      image = 'tmp/imgW-'+self.mtag+stag
     
      #self.getWImage(iexp,Timage,image)

      Rrwfl = eicop2d(swfl,gg,self.par) 
      Rrwfl.ADJT(adjR,res)

      #Rrwflim = cicop2d(swfl,self.par)
      #Rrwflim.ADJT(adjRim,image)

      Rswfl = eicop2d(rwfl,gg,self.par,positive='n') 
      Rswfl.ADJT(adjS,res)

      #Rswflim = cicop2d(rwfl,self.par)
      #Rswflim.ADJT(adjSim,image)

      L = aweop2d(vel,self.par,'','','',
                  custom=self.customfd)

      L.ADJT(aswfl,adjS)
      #L.ADJT(asimp,adjSim)

      L.FORW(adjR,arwfl)
      #L.FORW(adjRim,arimp)

      #self.corrgrad(sgradIm,swfl,asimp)
      #sgradsIm.append(sgradIm)

      #self.corrgrad(rgradIm,rwfl,arimp)
      #rgradsIm.append(rgradIm)

      self.corrgrad(sgrad,swfl,aswfl)
      sgrads.append(sgrad)

      self.corrgrad(rgrad,rwfl,arwfl)
      rgrads.append(rgrad)

      Iterate()
    Join()

    self.spartialgrad= 'tmp/spartial-grad-'+self.mtag
    self.rpartialgrad= 'tmp/rpartial-grad-'+self.mtag
    self.partialgrad= 'images/partial-grad-'+self.mtag

    #self.rpartialgradIM= 'tmp/rpartialIm-grad-'+self.mtag
    #self.spartialgradIM= 'tmp/spartialIm-grad-'+self.mtag
    #self.partialgradIm= 'tmp/partialIm-grad-'+self.mtag

    patternS = 'tmp/sgrad-'+self.mtag+'-%03d'
    patternR = 'tmp/rgrad-'+self.mtag+'-%03d'
    #patternSIm = 'tmp/sgradIm-'+self.mtag+'-%03d'
    #patternRIm = 'tmp/rgradIm-'+self.mtag+'-%03d'

    #self.mpiStack(self.spartialgradIM,sgradIm,pattern=patternSIm)
    #self.mpiStack(self.rpartialgradIM,rgradIm,pattern=patternRIm)

    self.mpiStack(self.spartialgrad,sgrad,pattern=patternS)
    self.mpiStack(self.rpartialgrad,rgrad,pattern=patternR)

    Flow(self.partialgrad,[self.spartialgrad,self.rpartialgrad],
      'add ${SOURCES[1:2]}')

    #Flow(self.partialgradIm,[self.spartialgradIM,self.rpartialgradIM],
    #  'add ${SOURCES[1:2]}')



  def mpiStack(self,image,images,pattern='tmp/img-%03d'):
    mpi.gridandstack(image,images,
     8,pattern,
     self.par['nz'],self.par['oz'],self.par['dz'],
     self.par['nx'],self.par['ox'],self.par['dx'],
     1.0         ,0.0         ,0.00225,
     nf=self.par['ns'],of=0,jf=1)

  def setScalar(self,scalar=1):
    self.scalar=scalar


  def getVel(self,isou,ivel,ovel,lower=1.45,upper=4.06):


    Flow(self.partialgrad,sgrads+rgrads,'add ${SOURCES[1:-1]} ')
    Flow(self.partialgradIm,sgradsIm+rgradsIm,'add ${SOURCES[1:-1]} ')

  def gradient(self,GR,OF,m):
    self.tag(m)
    self.of(OF,m)
    self.adjointStateVariables(m)
    Flow(GR,[self.partialgrad,m,self.mask],
    'math m=${SOURCES[1]} output="%g*(-2)/m^3*x1*(input)"|'%self.scalar+
     self.prec+'|add mode=p ${SOURCES[2]}')

    #Flow(GR,[self.partialgrad,self.partialgradIm,m],
    #'math m=${SOURCES[2]} im=${SOURCES[1]} output="x1*2/m^3*(input-%g*im)"|scale axis=123|'%self.im+self.prec)

  def corrgrad(self,grad,wfl1,wfl2):
    addition = ''
    if not self.minimize: addition = '|scale rscale=-1'
    Flow(grad,[wfl1,wfl2],
      '''
      nderiv order=2 axis=3 length=11 |
      add ${SOURCES[1]} mode=p|
      stack axis=3 norm=n 
      '''+addition)

  def getLocations(self,gg,image):
    imgtag = image.replace("/","_")
    msk = 'tmp/'+ imgtag
    nimage = 'tmp/'+ imgtag+'noise'
    Flow(msk,self.parcip['exclusion'],'costaper nw1=30 nw2=60')
    Flow(nimage,image,'scale axis=123 |scale rscale=100|noise | agc rect1=80')
    ciputils.pickcip(gg,nimage,msk,self.parcip)


  def tag(self,m):
    self.mtag = m.replace("/","_")

  def getVel(self,isou,ivel,ovel,lower=1.45,upper=4.06):
    Flow(ovel,ivel,
      '''
      remap1 order=1 n1=%(nz)d d1=%(dz)g o1=%(oz)g |'''%self.ispar(isou)+'''
      clip2 lower=%g upper=%g|
      transp |
      '''%(lower,upper)+''' 
      remap1 order=1 n1=%(nx)d d1=%(dx)g o1=%(ox)g |'''%self.ispar(isou)+
      '''
      clip2 lower=%g upper=%g|
      transp
      '''%(lower,upper))

  def getWImage(self,isou,ivel,ovel):
    Flow(ovel,ivel,
      '''
      remap1 order=1 n1=%(nz)d d1=%(dz)g o1=%(oz)g |
      transp |
      remap1 order=1 n1=%(nx)d d1=%(dx)g o1=%(ox)g |
      transp'''%self.ispar(isou))



  def ispar(self,isou):
    par = self.par
    sx = (isou)*par['ds'] + par['os'] 
    maxsx = int((sx + (par['maxh']+par['aper'])-par['ox'])/par['dx'])*par['dx']
    minsx = int((sx - par['aper']-par['ox'])/par['dx'])*par['dx']
    minsx = max(minsx,par['ox'])
    maxsx = min(maxsx,par['mx'])
    nx = int((maxsx-minsx)/par['dx'])
    opar = {'ox':minsx,'dx':par['dx'],'nx':nx,
            'oz':par['oz'],'dz':par['dz'],'nz':par['nz'],
            'nb':par['nb'],'jsnap':par['jsnap']}
    fdmod.param(opar)
    return opar


class ApertureRtm:
  def __init__(self,dS,ss,dR,rr,par):
    self.par = par.copy()
 
    self.dS = dS
    self.ss = ss
    self.nexp = len(ss)
    self.dR = dR
    self.rr = rr
    self.setCompute()
    self.customfd = 'order=25 jsnap=1 nb=20'

  def setCompute(self,nodes=10,time=5):
    self.nodes = nodes
    self.ipn =  int(self.nexp/nodes)
    self.time = time
    self.cluster = True


  def Rtm(self,ivel,image):
    images=[]
    vtag = ivel.replace("/","-")
    Fork(time=self.time,ipn=self.ipn,nodes=self.nodes)
    for isou in range(self.nexp):
      stag = '-'+vtag+'-%03d'%isou
      wav = self.dS[isou]
      ss  = self.ss[isou]
      dat = self.dR[isou]
      rr  = self.rr[isou]

      # outputs 
      swfl = 'tmp/swfl'+stag
      rwfl = 'tmp/rwfl'+stag
      img  = 'tmp/img'+stag
      vel  = 'tmp/vel'+stag

      # restrict domain:
      self.getVel(isou,ivel,vel)

      # operators
      Lswfl = aweop2d(vel,self.par,ss,'','',custom=self.customfd)
      Lrwfl = aweop2d(vel,self.par,'',rr,'',custom=self.customfd)

      Lswfl.FORW(wav,swfl)
      Lrwfl.ADJT(rwfl,dat)
    
      Rop = cicop2d(swfl,self.par)
      Rop.FORW(rwfl,img)
      images.append(img)
      Iterate()
    Join()

    mpi.gridandstack(image,images,
     8,'tmp/img-'+vtag+'-%03d',
     self.par['nz'],self.par['oz'],self.par['dz'],
     self.par['nx'],self.par['ox'],self.par['dx'],
     1.0         ,0.0         , 0.1,
     nf=self.par['ns'],of=0,jf=1)


  def Ximage(self,ivel,oimage,oximage,gg):
    self.Rtm(ivel,oimage)
    vtag = ivel.replace("/","-")

    eimages=[]
    Fork(time=self.time,ipn=self.ipn,nodes=self.nodes)
    for isou in range(self.nexp):
      stag = '-'+vtag+'-%03d'%isou
      swfl = 'tmp/swfl'+stag
      rwfl = 'tmp/rwfl'+stag
      eimg  = 'tmp/ximg'+stag
      ERop = eicop2d(swfl,gg,self.par)
      ERop.FORW(rwfl,eimg)
      eimages.append(eimg)
      Iterate()
    Join()

    Flow(oximage,eimages,'add ${SOURCES[1:-1]}')




  def getVel(self,isou,ivel,ovel,lower=1.45,upper=4.06):
    Flow(ovel,ivel,
      '''
      remap1 order=1 n1=%(nz)d d1=%(dz)g o1=%(oz)g |'''%self.ispar(isou)+'''
      clip2 lower=%g upper=%g|
      transp |
      '''%(lower,upper)+''' 
      remap1 order=1 n1=%(nx)d d1=%(dx)g o1=%(ox)g |'''%self.ispar(isou)+
      '''
      clip2 lower=%g upper=%g|
      transp
      '''%(lower,upper))

  def ispar(self,isou):
    par = self.par
    sx = (isou)*par['ds'] + par['os'] 
    maxsx = int((sx + (par['maxh']+par['aper'])-par['ox'])/par['dx'])*par['dx']
    minsx = int((sx - par['aper']-par['ox'])/par['dx'])*par['dx']
  
    minsx = max(minsx,par['ox'])
    maxsx = min(maxsx,par['mx'])
    nx = int((maxsx-minsx)/par['dx'])
    opar = {'ox':minsx,'dx':par['dx'],'nx':nx,
            'oz':par['oz'],'dz':par['dz'],'nz':par['nz'],
            'nb':par['nb'],'jsnap':par['jsnap']}
    fdmod.param(opar)
    return opar


  def getLocations(self,gg,image):
    imgtag = image.replace("/","_")
    msk = 'tmp/'+ imgtag
    nimage = 'tmp/'+ imgtag+'noise'
    Flow(msk,self.parcip['exclusion'],'costaper nw1=30 nw2=60')
    Flow(nimage,image,'scale axis=123 |scale rscale=100|noise | agc rect1=80')
    ciputils.pickcip(gg,nimage,msk,self.parcip)

  def setPicker(self,parcip):
    self.cip = True
    self.parcip = parcip.copy()


class wetIDMax(wetID):
  '''
    min -||H*Eimage|| 
  '''
  def setSigmas(self,sz,sx,st):
    self.st = st
    self.sx = sx
    self.sz = sz

  def computePenalty(self,m):
    # over-writes method
    gauss = 'math output="exp(-(x1^2/(2*%g^2)))*exp(-(x2^2/(2*%g^2)))*exp(-(x3^2/(2*%g^2))) " '%(self.sz,self.sx,self.st)

    self.penalty='tmp/penalty'+self.mtag
      
    Flow(self.penalty,[m,self.gg],
      '''
      python ./getpenalty.py nht=%(nht)d nhx=%(nhx)d nhz=%(nhz)d 
                             dht=%(dt)g  dhx=%(dx)g  dhz=%(dz)g
             cip=${SOURCES[1]}|
      '''%self.par+
      gauss)

    self.residual = 'images/res-'+self.mtag
    self.res = 'images/resx-'+self.mtag
    sources = [self.eimage,self.penalty]

    Flow([self.residual,self.res],sources,
      'add mode=p ${SOURCES[1]}> ${TARGETS[1]} && add mode=p ${SOURCES[1]} <${TARGETS[1]}|scale rscale=-1')
  
  def of(self,OF,m):
    self.tag(m)

    self.stateVariables(m)
    self.computePenalty(m)

    Flow(OF,[self.res],
       'add mode=p $SOURCE |stack axis=0 norm=n |scale rscale=-%g'%(self.scalar/self.parcip['nmax']))

  def gradient(self,GR,OF,m):
    self.tag(m)
    self.of(OF,m)
    self.adjointStateVariables(m)
    Flow(GR,[self.partialgrad,m,self.mask],
    'math m=${SOURCES[1]} output="%g*(2)/m^3*x1^2*(input)"|'%self.scalar+
     self.prec+'|add mode=p ${SOURCES[2]}')
