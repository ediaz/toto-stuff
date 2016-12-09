try: 
  from rsf.cluster import *
  from rsf.proj import WhereIs
except: from rsf.proj import *

from aweOP import aweop2d
from injOP import injop2d 

java = WhereIs('java')+' -server -ea -Xmx12g '
ediaz = ' ediaz.rsf.'
lcorr = java+ediaz+'LocalCorr '



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
      lcorr+' g=${SOURCES[1]} s3=%g nlag=%d fwd=y'%(self.s3,self.nlag)
      )

  def ADJT(self,f,corr):
    Flow(f,[corr,self.g],
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

  def setPrec(self,rect1=5,rect2=5,repeat=5):
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
      mdatT = 'corrs/mdat-'+self.mtag+stag
      mdat = 'tmp/mdat-'+self.mtag+stag
      dat  = 'tmp/dat-'+self.mtag+stag
      corr = self.sdir+'corr-'+self.mtag+stag

      Ldat = aweop2d(vel,self.par,ss,rr,'',custom=self.customfd)
      Ldat.FORW(wav,mdatT)

      if not feval: 
        Lswfl = aweop2d(vel,self.par,ss,'','',
                        custom=self.customfd)
        Lswfl.FORW(wav,swfl)


      Flow(mdat,mdatT,'costaper nw1=1 |transp')
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

      Flow(pcorr,[corr,self.penalty],
        '''
        python applyPenalty.py pen=${SOURCES[1]} 
        ''')
      Flow(ofn,[corr,self.penalty],
        '''
        python applyPenalty.py pen=${SOURCES[1]} |
        math output="input^2"|
        stack axis=0 norm=n''')


      Iterate()
    Join()

    ofn = 'tmp/ofn-'+self.mtag
    Flow(OF,map(lambda x:'tmp/ofn-'+self.mtag+'-%03d'%x,range(self.nexp)),
      'cat ${SOURCES[1:%d]}|stack axis=0 norm=n '%self.nexp)
    self.OF  = OF
    


  def adjointStateVariables(self,m,feval=False):
    self.tag(m)
    zero = 'tmp/zero'+self.mtag
    Flow(zero,'tmp/dat-'+self.mtag+'-000','math output="0"')
    
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


      nadjtsrcT = 'tmp/nadjtsrcT-'+self.mtag+stag
      dadjtsrcT = 'tmp/dadjtsrcT-'+self.mtag+stag
      adjsrc = 'tmp/adjtsrc-'+self.mtag+stag
      
      grds = 'tmp/grd-'+self.mtag+stag
  
      lcorr = Lcorr(dat,self.nlag,self.s3,
                    self.par['ot'],self.par['dt'],self.par['nt'],
                    self.par['ot'],self.par['dt'],self.par['nt'])

      lcorr.ADJT(nadjtsrcT,pcorr)
      

      Flow(adjsrc,[nadjtsrcT],'transp')

      Lrwfl = aweop2d(vel,self.par,'',rr,'',
                        custom=self.customfd)
      Lrwfl.ADJT(awfl,adjsrc)

      Flow(grds,[swfl,awfl],
        'nderiv length=3 axis=3 order=2 |add mode=p ${SOURCES[1]} |stack axis=3 ')

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
    'math m=${SOURCES[1]} output="-2e7/m^3*input"|'+self.prec+'|scale axis=123')




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



