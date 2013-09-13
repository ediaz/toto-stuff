#!/usr/bin/env python

import stopwatch as sw
import numpy as np 
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import lab4
 




def plot(img,fig=1):
  plt.figure(fig)
  fig = plt.imshow(img,cmap=cm.gist_gray)

  
def getimage():
  tmp=mpimg.imread('esteban.png')

  img = like(tmp)
  img = tmp

  return img

def getdim(img):
  n2 = len(img[0])
  n1 = len(img) 
  return n1,n2

def like(x):
  n1,n2 = getdim(x)
  return np.zeros((n1,n2),dtype='float32',order='F')

 
def goSmooth1():
  x = getimage() 
  a = float(0.9)
  y = like(x)

  n1,n2 = getdim(x) 
  plot(x,1)

  lab4.dsp.smooth1(a,x,y,n1,n2)
  print "   mean = %13.11f  "%lab4.dsp.mean(y,n1*n2)
  plot(y,2)

def goSmooth2():
  x = getimage() 
  a = float(0.9)
  y = like(x)

  n1,n2 = getdim(x) 

  lab4.dsp.smooth2(a,x,y,n1,n2)
  plot(y,3)
  
  print "   mean = %13.11f "%lab4.dsp.mean(y,n1*n2)

def goSmooth2T():
  x = getimage() 
  a = float(0.9)
  y = like(x)

  n1,n2 = getdim(x) 

  lab4.dsp.smooth2t(a,x,y,n1,n2)
  plot(y,5)
  
  print "   mean = %13.11f "%lab4.dsp.mean(y,n1*n2)



def goSmooth():
  x = getimage() 
  a = float(0.9)
  y = like(x)

  n1,n2 = getdim(x) 

  lab4.dsp.smooth1(a,x,y,n1,n2)
  lab4.dsp.smooth2(a,y,y,n1,n2)
  plot(y,4)
  
  print "   mean = %13.11f "%lab4.dsp.mean(y,n1*n2)
  

goSmooth1()
goSmooth2()
goSmooth()


'''
print "Python-Fortran"

n = 1001
maxtime = 2.0 
a =float(0.99) 

x= np.zeros(n, dtype='float32')

x[0]=x[int(n/2)]= x[n-1]= 1.0 

nsmooth=0
y= np.zeros(n, dtype='float32')

# initialize stopwatch object:
t= sw.stopwatch()


t.start()
while t.time() < maxtime :
    lab4.dsp.smooth_dsp(a,x,y)
    nsmooth+=1
t.stop()

print "nsmooth = %d "%nsmooth
print "   mean = %9.7f "%lab4.dsp.mean_dsp(y)
print "   time = %9.7f "%t.time()
print " mflops = %f "%((6.0e-6*n*nsmooth/t.time()))
'''


plt.show()
