#!/usr/bin/env python

import gaussmix as gm
import numpy as np
from scipy import stats
import steamleaf as stl


n=30


'''
#plot i
x1= gm.mix_2gauss(n,0,0,1,36,.95)

gm.histogram_plot(x1,200,'fig1.png')


#plot ii
x2= gm.mix_2gauss(n,0,3,1,1,.80)
#gm.histogram_plot(x2,80,'fig2.png')

#plot iii
x3= gm.lognormal(0,1,n)
#gm.histogram_plot(x3,100,'fig3.png')
'''



# Problem 4.11
data = np.array( [ 97.0,97.2,97.3,97.6,97.6,97.7,97.9,98.2,98.2,98.4,
                   98.4,98.5,98.6,98.6,98.6,98.6,98.6,98.7,98.8,98.9,
                   99.0,99.0,99.1,99.2,99.3,99.4,99.5,99.5,99.7,99.8])


print 'Problem 4: mean   = ',data.mean()
print 'Problem 4: median = ',np.median(data)
print 'Problem 4: std    = ',data.std()
print 'Problem 4: P25,P50,P75    = ',stats.scoreatpercentile(data, 25),stats.scoreatpercentile(data, 50),stats.scoreatpercentile(data, 75)



print '=================='
print '=  stem and leaf ='
print '=================='

sl = stl.stemleafpairs(data, stempos = 0, leafwidth = 1)
stl.prettyprint(sl,stemwidth = 3)


# Problem 4.13

rain = np.array( [1468,  909,  841,  475,  846,  452,
                  3830, 1397,  556,  978, 1715,  747,
                   909, 2002, 1331, 1227, 2543, 2649,
                  1781, 1717, 2718,  584, 1859, 1138,
                  2675, 1872, 1359, 1544, 1372, 1334,
                   955, 1849,  719, 1737, 1389,  681,
                  1565,  701,  994, 1188,  962, 1564,
                  1800,  580, 1106,  880,  850] )

print '=================='
print '=  stem and leaf ='
print '=================='

sl = stl.stemleafpairs((rain)*0.005, stempos = 0 , leafwidth = 1)
stl.prettyprint(sl,stemwidth = 3)



#Box plot

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10,6))
fig.canvas.set_window_title('A Boxplot Example')
bp = plt.boxplot(rain, notch=0, sym='+', vert=1, whis=1.5)
ax1 = fig.add_subplot(111)
# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.5)

plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='red', marker='+')


# QQ-plot :

import qqplot as qq

qq.qqplots(rain)


plt.show()
#plot iii
