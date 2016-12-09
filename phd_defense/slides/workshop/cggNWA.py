try:
    from rsf.cluster import *
except:
    from rsf.proj import *
from rsf.recipes import fdmod
import os
# ------------------------------------------------------------
# This module reads the data, extract geometry, creates
# the source file (provided by CGG)
#
# You can extract the navegation data with the function 
# getxy()
#
# If you do that you can see that the last shots gather
# streamer have feathering. 
#
# The function project takes the original xy data and 
# 'projects' it into a line, which the user can define using (xo,yo,angle)
#
# I assumed that the data was acquired in a perfect 2D sense,
# if you do that you will get a conventional RSF cube
# axis 1: time
# axis 2: offset
# axis 3: shot position
#
# This is the easiest way to bin the data, it might not be
# the best one...
#
# Esteban Diaz
# Center for Wave Phenomena
# 2013.9
# ------------------------------------------------------------

ddir = '%s/data/NWA006dist/'%os.getenv('CWP')
# ------------------------------------------------------------
# model parameters
def getpar():
    par = {
        'nx':6119, 'ox':0.0000, 'dx':0.00625,  'lx':'x',      'ux':'km',
        'nz':2401,  'oz':0.0000, 'dz':0.00625,  'lz':'z',      'uz':'km',
        'nt':3526,  'ot':0.0000, 'dt':0.00200,  'lt':'t',      'ut':'s',
        'nh':648,   'oh':0.1690, 'dh':0.01250,  'lh':'offset', 'uh':'km',
        'ns':1824,  'os':0.0   , 'ds':0.01875,  'ls':'sx',     'us':'km',
        'sdepth':0.005,
        'rdepth':0.0125
        }
    par['nb']=250
    return par

def getfile(which='broadband'):
  files = {'broadband':'NW_Australia_BroadSeis_SP_Line006_CGGConfidential.segy',
           'conventional':'NW_Australia_Conventional_SP_Line006_CGGConfidential.segy'}
  fname = files[which] 
  ftag = fname.replace("_SP_Line006_CGGConfidential.segy","")
  return ftag,fname 

def SegyShots2rsf(par,which='broadband'):
  ftag,fname = getfile(which)

  header = 'headers_'+ftag
  ebcdic = 'EBCDIC_'+ftag
  
  Flow(ftag+'.ebcdic',ddir+fname,'\dd if=$SOURCE count=40 bs=80 conv=unblock,ascii',stdin=0)
  Flow([ftag,header],ddir+fname,
    '''
    segyread tfile=${TARGETS[1]} |
    put n2=%(nh)d o2=%(oh)g d2=%(dh)g 
        n3=%(ns)d o3=%(os)g d3=%(ds)g                                
    '''%par)  

def extractShot(par,which='broadband',snum=785,stag='shot-',sstag='ss-',rrtag='rr-'):
  ftag,fname = getfile(which)
  Flow(stag+'%04d'%snum,[ftag],
    '''
    window n3=1 f3=%d
    '''%snum) 

def getchannel(chansection,head,data,chan):  
  # this function computes a common channel section, 
  # useful for extracting near, mid, far channel
  Flow(chansection,[data],
    '''
    window n2=1 f2=%d
    '''%(chan-1))

def getxy(out,head,which='source'):  
  s = {'source':21,'receiver':23}
  f1 = s[which]
  par = getpar() 
  Flow(out,head,'window n1=2 f1=%d |dd type=float | scale rscale=0.001 |'%f1+'put n2=%(nh)d n3=%(ns)d'%par)

def plot_globalxy(col=6):
  return '''
         window n1=2 | 
         dd type=complex |
         graph labelrot=n wantaxis=y title=" "  symbol=. plotcol=%d plotfat=10 
               label1="x" unit1="km" label2="y" unit2="km" 
               min1=265 max1=310 
               min2=7755 max2=7775 
               screenratio=%g screenht=5 yll=1
         '''%(col,(7775.-7755)/(310.-265.))


def project(out,ifile,alpha,xo,yo):
  from math import pi 
  '''
  alpha: azimuth angle of the new x axis
  xo, yo : coordinates of the new transformation 0
  '''
  xtmp = 'xtmp_'+ifile 
  ytmp = 'ytmp_'+ifile 

  xotmp = 'xotmp_'+ifile 
  yotmp = 'yotmp_'+ifile 

  Flow(out,ifile,
    'window n1=1 squeeze=n > %s; '%(xtmp)+ 
    '$RSFROOT/bin/sfwindow <$SOURCE n1=1 f1=1 squeeze=n >%s;'%ytmp+
    '$RSFROOT/bin/sfmath x=%s y=%s output="(x-%g)*cos(%g) -(y-%g)*sin(%g) " > %s;'%(xtmp,ytmp,xo,alpha*pi/180,
                                                                                   yo,alpha*pi/180,xotmp)+
    '$RSFROOT/bin/sfmath x=%s y=%s output="(x-%g)*sin(%g) +(y-%g)*cos(%g) " > %s;'%(xtmp,ytmp,xo,alpha*pi/180,
                                                                                   yo,alpha*pi/180,yotmp)+
     '$RSFROOT/bin/sfcat %s %s axis=1 >$TARGET;'%(xotmp,yotmp)+
     '$RSFROOT/bin/sfrm %s %s %s %s ; echo '%(xotmp,yotmp,xtmp,ytmp)
     ,stdout=0)

def getcable(cable,header,par):
  '''
  This function computes the average streamer cable
  '''
  Flow(cable,header,
    '''
    window n1=1 f1=18 squeeze=n|dd type=float|
    scale rscale=0.0001 |
    put n2=%(nh)d n3=%(ns)d |
    window n3=300 squeeze=n |
    stack axis=3 |
    window |
    put o1=%(oh)g d1=%(dh)g n1=%(nh)d |
    rtoc|
    math output="x1+I*input"|transp|
    dd type=float 
    '''%par)


def getSource(wav,par):
  source = ddir+'04*.txt'
  Flow(wav,None,
    '''
    \cat %s > tmp$TARGET;
    echo n1=501 d1=0.002 o1=0 >> tmp$TARGET; 
    echo data_format="ascii_float" >>tmp$TARGET;
    echo in=tmp$TARGET >>tmp$TARGET;
    <tmp$TARGET  $RSFROOT/bin/sfdd form=native |pad n1=4000|window n1=%d >$TARGET ;
    \\rm tmp$TARGET
    '''%(source,par['nt']),stdin=0,stdout=2)
