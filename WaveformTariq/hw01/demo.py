try:    from rsf.cluster import *
except: from rsf.proj    import *
from rsf.recipes import wplot,awe,geom

def param():
    par = dict(nt=2500, ot=-1.000, dt=0.002, lt='t', ut='s',  nht=500,
               nx=301,  ox=-1.000, dx=0.020, lx='x', ux='km', nhx=20,
               ny=1,    oy=0,      dy=0.000, ly='y', uy='km', nhy=0,
               nz=181,  oz=-0.200, dz=0.020, lz='z', uz='km', nhz=20,
               kt=500,nb=100,frq=10,jsnap=50,
               gaus='y',verb='y',
               vel=2.0,vsc=1.02)

    par['shft']=35 # shift source/receivers from the edge
    par['oz']=-par['shft']*par['dz']

    par['zrefl']=2.0

    par['tfrm']=3.25 # max movie time
    par['nfrm']=9    # number of frames
    par['jfrm']=int((par['tfrm']/par['dt'])/par['nfrm'])
    par['ffrm']=par['kt']

    par['ixcip']=(4.0-par['ox'])/par['dx']; xcip=par['ox']+par['ixcip']*par['dx']
    par['izcip']=par['shft'];               zcip=par['oz']+par['izcip']*par['dz']
    
    return par

# ------------------------------------------------------------
def myframes(file,byte,custom,par):
    Flow('S-'+file, file,'window')
    wplot.inine('S-'+file,byte,custom,par,scale=0.4,ymax=6.5,ratio=0.65)
    
# ------------------------------------------------------------
def wavelet(wav,par):
    awe.wavelet(wav,par['frq'],'',par)
    Result(wav,'transp |' + wplot.waveplot('',par))

# ------------------------------------------------------------
# point source
def ss(ss,par):
    par['ixsou']=(0.0-par['ox'])/par['dx']; xsou=par['ox']+par['ixsou']*par['dx']
    par['izsou']=par['shft'];               zsou=par['oz']+par['izsou']*par['dz']
    geom.point2d(ss,xsou,zsou,'',par)
    Plot(ss,wplot.rrplot2d('plotcol=1 plotfat=20',par))

# ------------------------------------------------------------
# point receiver
def rr(rr,par):
    par['ixrec']=(4.0-par['ox'])/par['dx']; xrec=par['ox']+par['ixrec']*par['dx']
    par['izrec']=par['shft'];               zrec=par['oz']+par['izrec']*par['dz']
    geom.point2d(rr,xrec,zrec,'',par)
    Plot(rr,wplot.rrplot2d('plotcol=2 plotfat=20',par))
    
# ------------------------------------------------------------
def vtran(vel,par):
    geom.zero2d(vel+'zero','',par)
    Flow(  vel,vel+'zero','math output="%(vel)g"'%par)
    Plot(  vel, wplot.igrey2d('mean=y',par))

def vback(vbk,vel,scale,par):
    Flow(vbk,vel,'scale rscale=%g'%scale)
    Plot(vbk,wplot.igrey2d('mean=y',par))
    
def vrefl(vel,par):
    geom.zero2d(vel+'zero','',par)
    Flow(vel,vel+'zero',
         '''
         spike nsp=1 mag=-%g k1=%d l1=%d |
         add add=%g
         '''%(par['vel'],(par['zrefl']-par['oz'])/par['dz']+1,par['nz'],par['vel']))
    Plot(vel,wplot.igrey2d('mean=y',par))


