try: 
  from rsf.cluster import *
  from rsf.proj import WhereIs
except: from rsf.proj import *

from aweOP import aweop2d
from injOP import injop2d 

java = WhereIs('java')+' -server -ea -Xmx12g '
ediaz = ' ediaz.rsf.'
lcorr = java+ediaz+'LocalCorr '



import os

AWE2d = os.getenv('iTeam')+'/m8rSolver/CODE/AWEOP2D.x'
AWE3d = os.getenv('iTeam')+'/m8rSolver/CODE/AWEOP3D.x'
#AWE2d = 'aweop2d'
#AWE3d = 'aweop3d'

# ------------------------------------------------------------
def awepar(par):
    if(not par.has_key('dabc')):     par['dabc']='n'
    if(not par.has_key('nb')):       par['nb']=0
    if(not par.has_key('ompchunk')): par['ompchunk']=1
    if(not par.has_key('ompnth')):   par['ompnth']=0
    if(not par.has_key('fsrf')):     par['fsrf']='n'
    if(not par.has_key('verb')):     par['verb']='n'
    awe = ' ' + \
        '''
            ompchunk=%(ompchunk)d ompnth=%(ompnth)d
            verb=%(verb)s fsrf=%(fsrf)s
            dabc=%(dabc)s nb=%(nb)d
            '''%par + ' '
    return awe

# ------------------------------------------------------------
class aweop2d:
    def __init__(self,vel,par,ss,rr,den,custom=''):
        self.vel=vel
        self.par=par.copy() # copy the par file instead of referencing it
        self.custom=custom
        
        awepar(self.par) # set AWE defaults
    
        self.dep = [self.vel]
    
        self.ss=''
        if(ss!=''):
            self.ss=' sou='+ss+'.rsf '
            self.dep.append(ss)

        self.rr=''
        if(rr!=''):
            self.rr=' rec='+rr+'.rsf '
            self.dep.append(rr)

        self.den=''
        if(den!=''):
            self.den=' den='+den+'.rsf '
            self.dep.append(den)
    
    def FORW(self,m,d):
        Flow(d,[m]+self.dep,
             AWE2d + ' adj=n vel=${SOURCES[1]}'+self.ss+self.rr+self.den
             + awepar(self.par)
             +' '+self.custom+
             '''
             > ${TARGET}_tmp &&
             transp <${TARGET}_tmp |
             halfint inv=y |transp >$TARGET &&
             rm ${TARGET}_tmp 
             ''',stdout=-1)
    
    def ADJT(self,m,d):
        Flow(m,[d]+self.dep,
             '''
             transp |
             halfint inv=y adj=y|transp >${TARGET}_tmp&&             
             < ${TARGET}_tmp
             '''+
             AWE2d + ' adj=y vel=${SOURCES[1]}'+self.ss+self.rr+self.den
             + awepar(self.par)
             +' '+self.custom+'>$TARGET&& rm ${TARGET}_tmp',stdout=-1)


class Lcorr:
  '''
  Local correlation object
  '''
  def __init__(self,g,nlag,s3,
                    ot,dt,nt,
                    otc,dtc,ntc):
    self.g = g
    self.nlag = nlag
    self.s3 = s3
    self.window = 'sinc o1=%g d1=%g n1=%d '%(otc,dtc,ntc)
    self.interp = 'sinc o1=%g d1=%g n1=%d '%(ot,dt,nt)

  def FORW(self,f,corr):
    Flow(corr,[f,self.g],
      lcorr+' g=${SOURCES[1]} s3=%g nlag=%d fwd=y|'%(self.s3,self.nlag)+
      self.window)

  def ADJT(self,f,corr):
    Flow(f,[corr,self.g],
      self.interp+'|'+
      lcorr+' g=${SOURCES[1]} s3=%g nlag=%d fwd=n'%(self.s3,self.nlag))







class LocalNormCorr:
  '''
    This object optimizes:
    
    J = || T Lc g || / || Lc g|| 
  '''
  def __init__(self,dS,ss,dR,rr,nlag,s3,par,penalty,customfd='jsnap=10 nb=100'):
    self.par = par.copy()
    
    self.dS = dS
    self.ss = ss
  
    self.dR = dR 
    self.rr = rr

    self.nexp = len(ss)

    self.nlag = nlag
    self.s3  = s3

    self.penalty = penalty

    self.customfd = customfd
    self.setSaveDir()
    self.setPrec()
    self.setCompute()

  def setPrec(self,rect1=1,rect2=1,repeat=3):
    self.prec = 'smooth rect1=%d rect2=%d repeat=%d'%(rect1,rect2,repeat)

  def setCompute(self,nodes=1):
    self.nodes = nodes
    self.ipn = int(self.nexp/nodes)

  def setSaveDir(self,sdir='corrs/'):
    self.sdir = sdir


  def stateVariables(self,m,feval=False):
    self.tag(m)

    Fork(time=4,ipn=self.ipn,nodes=self.nodes)
    for iexp in range(self.nexp):
      stag = '-%03d'%iexp
      ss   = self.ss[iexp] 
      rr   = self.rr[iexp]
      wav  = self.dS[iexp]
      datT  = self.dR[iexp]
      vel  = m 
      swfl = 'tmp/swfl-'+self.mtag+stag
      mdatT = 'tmp/mdatT-'+self.mtag+stag
      mdat = 'tmp/mdat-'+self.mtag+stag
      dat  = 'tmp/dat-'+self.mtag+stag
      corr = self.sdir+'corr-'+self.mtag+stag

      Ldat = aweop2d(vel,self.par,ss,rr,'',custom=self.customfd)
      Ldat.FORW(wav,mdatT)

      if not feval: 
        Lswfl = aweop2d(vel,self.par,ss,'','',
                        custom=self.customfd)
        Lswfl.FORW(wav,swfl)


      Flow(mdat,mdatT,'transp')
      Flow(dat ,datT,'transp')
      lcorr = Lcorr(dat,self.nlag,self.s3,
                    self.par['ot'],self.par['dt'],self.par['nt'],
                    self.par['ot'],self.par['dt'],int(self.par['nt']))
      lcorr.FORW(mdat,corr)

      Iterate()
    Join()


  def of(self,OF,m,feval=False):
    self.stateVariables(m,feval)

    Fork(time=4,ipn=self.ipn,nodes=self.nodes)
    for iexp in range(self.nexp):
      stag = '-%03d'%iexp
      mdat = 'tmp/mdat-'+self.mtag+stag
      dat  = 'tmp/dat-'+self.mtag+stag
      corr = self.sdir+'corr-'+self.mtag+stag
      pcorr = 'tmp/pcorr-'+self.mtag+stag 
      ofn  = 'tmp/ofn-'+self.mtag+stag  # numerator norm
      ofd  = 'tmp/ofd-'+self.mtag+stag  # denomitator norm
      of  = 'tmp/of-'+self.mtag+stag  # denomitator norm

      Flow(pcorr,[corr,self.penalty],
        '''
        python applyPenalty.py pen=${SOURCES[1]} |
        python applyPenalty.py pen=${SOURCES[1]} 
        ''')
      Flow(ofn,[corr,self.penalty],
        '''
        python applyPenalty.py pen=${SOURCES[1]} |
        math output="input^2"|
        stack axis=0 norm=n''')

      Flow(ofd,corr,'math output="input^2" | stack axis=0 norm=n')
      Flow(of,[ofn,ofd],'add mode=d ${SOURCES[1]}')

      Iterate()
    Join()

    ofn = 'tmp/ofn-'+self.mtag
    Flow(ofn,map(lambda x:'tmp/ofn-'+self.mtag+'-%03d'%x,range(self.nexp)),
      'cat ${SOURCES[1:%d]}|stack axis=0 norm=n '%self.nexp)

    ofd = 'tmp/ofd-'+self.mtag
    Flow(ofd,map(lambda x:'tmp/ofd-'+self.mtag+'-%03d'%x,range(self.nexp)),
      'cat ${SOURCES[1:%d]}|stack axis=0 norm=n '%self.nexp)
    Flow(OF,[ofn,ofd],'add mode=d ${SOURCES[1]}')


  def adjointStateVariables(self,m,feval=False):
    self.tag(m)

    Fork(time=4,ipn=self.ipn,nodes=self.nodes)
    for iexp in range(self.nexp):
      stag = '-%03d'%iexp
      ss   = self.ss[iexp] 
      rr   = self.rr[iexp]
      wav  = self.dS[iexp]
      datT  = self.dR[iexp]
      vel  = m 
      swfl = 'tmp/swfl-'+self.mtag+stag
      awfl = 'tmp/awfl-'+self.mtag+stag
      dat  = 'tmp/dat-'+self.mtag+stag
      pcorr = 'tmp/pcorr-'+self.mtag+stag
      corr = self.sdir+'corr-'+self.mtag+stag

      of = 'tmp/of-'+self.mtag+stag  # denomitator norm
      nof = 'tmp/nof-'+self.mtag+stag
    
      ofn  = 'tmp/ofn-'+self.mtag+stag  # numerator norm
      nofd  = 'tmp/nofd-'+self.mtag+stag  # numerator norm
      ofd  = 'tmp/ofd-'+self.mtag+stag  # denomitator norm
      dnofd  = 'tmp/dnofd-'+self.mtag+stag  # denomitator norm
      zero = ofn+'-zero'

      nadjtsrcT = 'tmp/nadjtsrcT-'+self.mtag+stag
      dadjtsrcT = 'tmp/dadjtsrcT-'+self.mtag+stag
      adjsrc = 'tmp/adjtsrc-'+self.mtag+stag
      Flow(nof,of,'scale rscale=-1')
      Flow(dnofd,ofd,'math output="-1/input"')
      
      
      grds = 'tmp/grd-'+self.mtag+stag
  
      lcorr = Lcorr(dat,self.nlag,self.s3,
                    self.par['ot'],self.par['dt'],self.par['nt'],
                    self.par['ot'],self.par['dt'],1+int(self.par['nt']))

      lcorr.ADJT(nadjtsrcT,pcorr)
      lcorr.ADJT(dadjtsrcT,corr)
      
      Flow(zero,dadjtsrcT,'math output="0"')

      #Flow(adjsrc,[dadjtsrcT,nadjtsrcT,nof,dnofd,zero],
      #  '''
      #  axplusy afile=${SOURCES[2]} y=${SOURCES[1]}  | 
      #  axplusy afile=${SOURCES[3]} y=${SOURCES[4]}  |
      #  transp |
      #  math output="input" 
      #  ''')

      Flow(adjsrc,[nadjtsrcT],
        '''
        transp 
        ''')

      Lrwfl = aweop2d(vel,self.par,'',rr,'',
                        custom=self.customfd)
      Lrwfl.ADJT(awfl,adjsrc)

      Flow(grds,[swfl,awfl],
        'nderiv axis=3 order=2 |add mode=p ${SOURCES[1]} |stack axis=3 ')

      Iterate()
    Join()
    Flow('tmp/partial-grad'+self.mtag,
          map(lambda x:'tmp/grd-'+self.mtag+'-%03d'%x,range(self.nexp)),
          'add ${SOURCES[1:-1]}')


  def gradient(self,GR,OF,m):
    self.tag(m)
    self.of(OF,m)
    self.adjointStateVariables(m)
    Flow(GR,['tmp/partial-grad'+self.mtag,m],
    'math m=${SOURCES[1]} output="2e8/m^3*input"|'+self.prec)




  def tag(self,m):
    self.mtag = m.replace("/","_")




class LocalCorr:
  def __inti__(self):
    pass


class GlobalCorr:
  def __init__(self):
    pass


class GlocalNcorr:
  def __init__(self):
    pass



