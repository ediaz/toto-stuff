from library import *


def go_fig(mode,eps2,fn=1):
  s = {'slow':'SV-waves','fast':'P-waves'}

  layer1 = layer(0.5, 3.0, 1.8, 0.2, 0.1,mode)
  layer2 = layer(1.0, 3.0, 1.5, eps2,0.3,mode)
  layers = [layer1,layer2]
  
  fig = plt.figure(fn*2)
  fig2 = plt.figure(fn*2+1)

  ax = fig.add_subplot(111)
  ax.set_title('rays for %s, $\epsilon=%3.2g$'%(s[mode],eps2))

  ax2 = fig2.add_subplot(111) 
  ax2.set_title('travel times for %s, $\epsilon=%3.2g$'%(s[mode],eps2))

  theta = np.linspace(-85.,85.,1000) 
  i = 0 
  dray = 10 # how many other rays you want to plot?
  for t in theta:
    (rpath,time)=ray(layers,t)
    x = rpath[:,0]
    z = rpath[:,1]
    offset = x[len(x)-1]
    if i%dray ==0 and abs(offset)<2.:
      ax.plot(x,z,'k')
      tmax = t
    if (abs(offset)<2.):
      ax2.plot(offset,time,'k.')
    i +=1
  print tmax
  ax.set_xlabel("x(km)")
  ax.set_ylabel("z(km)")
  ax.set_xlim(-2.,2.)
  ax.invert_yaxis() 
  
  ax2.set_xlabel("x(km)")
  ax2.set_ylabel("t(s)")
  ax2.set_xlim(-2.,2.)
  #ax2.set_ylim(.8,3.2)
  ax2.invert_yaxis() 

  fig.savefig("report/Fig/rays_%s_eps%d.png"%(mode,int(eps2*10)))
  fig2.savefig("report/Fig/tt_%s_eps%d.png"%(mode,int(eps2*10)))

ifn=1
for mode in ['fast','slow']:
  for eps2 in [0.2,-0.2]:
    go_fig(mode,eps2,ifn)
    ifn+=1



