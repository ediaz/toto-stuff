#!/usr/bin/env python

import gaussmix as gm

n=1000




#plot i
x1= gm.mix_2gauss(n,0,0,1,36,.95)

gm.histogram_plot(x1,200,'fig1.png')


#plot ii
x2= gm.mix_2gauss(n,0,3,1,1,.80)
#gm.histogram_plot(x2,80,'fig2.png')

#plot iii
x3= gm.lognormal(0,1,n)
#gm.histogram_plot(x3,100,'fig3.png')


#plot iii
