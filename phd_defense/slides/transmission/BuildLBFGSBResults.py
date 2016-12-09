 #!/usr/bin/env python
import rsf.api as rsf
import numpy as np
import sys,os

par = rsf.Par()

FmodelMov = par.string('modelmov')
Fof = par.string('of')
Fgrad = par.string('gradmov')


dbase = par.string("dbasepath")  # this is where we put the final results
inversion = par.string("inversionpath",'')  # this from where
                                         # we get the results

of = np.zeros(1,'f')
try:
  A = np.loadtxt(dbase+'iter.log')
  try:
    it = A[:,6]
  except:
    it = [A[6]]
except:
  sys.exit(0)

ofl = ''
xl = ''
gl = '' 
for f in it:
  ofl += ' '+inversion+'of-%03d.rsf'%(f-1)
  xl  += ' '+inversion+'xf-%03d.rsf'%(f-1)
  gl  += ' '+inversion+'grad-%03d.rsf'%(f-1)



i=0
  
def cat(list_in,file_out):
  rsfroot=os.environ['RSFROOT']+'/bin/'
  global i
  i+=1
  if not file_out==None:
    command = rsfroot+'sfcat '+list_in+' axis=4 |'+ \
              rsfroot+'sfwindow out=stdout > '+file_out
    print >> sys.stderr, command
    os.system(command)
  else: 
    print >> sys.stderr,'list %d was not built'%i
    

cat(ofl,Fof)
cat(xl,FmodelMov)
cat(gl,Fgrad)
