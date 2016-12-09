try:    from rsf.cluster import *
except: from rsf.proj    import *
import rsflof,rsflsf,cippicker

def pickcip(pck,img,msk,parcip):

    # linearity
    rsflof.linearity2D(img+'-lin',
                       img,
                       parcip['sigma1'],parcip['sigma2'],
                       parcip['lzclip'],parcip['jmem'])

    # semblance
    rsflsf.linearSemblance2D(img+'-smb',
                             img,
                             parcip['hw1'],parcip['hw2'],
                             parcip['sigma1'],parcip['sigma2'],
                             parcip['szclip'],parcip['jmem'])

    # envelope
    Flow(img+'-env',img,'agc | envelope | smooth rect1=5 rect2=5')

    # priority map
    Flow(img+'-map',
         [img+'-lin',img+'-smb',img+'-env',msk],
         '''
         math output="i*j*k*l"
         i=${SOURCES[0]} j=${SOURCES[1]} k=${SOURCES[2]} l=${SOURCES[3]}
         ''',stdin=0)
    
    # picking
    cippicker.pickCoord2D(img+'-pck',
                          img+'-map',
                          parcip['r1'],parcip['r2'],parcip['jmem'])

    Flow(pck,
         img+'-pck','window n1=%d | transp |reverse which=1 opt=i'%(parcip['nmax']))
