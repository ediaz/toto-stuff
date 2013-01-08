#!/usr/bin/env python


def smooth (a, x, y):
    n = len(x) 
    b = 1.0-a 
    yi = y[0] = x[0]

    for i in range(1,n-1,1):
        y[i] = yi = a*yi + b*x[i] 

    y[n-1] = yi = (a*yi +x[n-1])/(1.0+a) 
    for i in range(n-2,-1,-1):
        y[i] = yi = a*yi + b*y[i]

def smooth2 (a, x, y):
    n = len(x) 
    b = 1.0-a
    sx = b ; sy = a 
    yi = y[0] = sx*x[0]

    for i in range(1,n-1,1):
        y[i] = yi = a*yi + b*x[i] 

    sx /= 1.0+a
    sy /= 1.0+a
    y[n-1] = yi = sy*yi+sx*x[n-1];

#    y[n-1] = yi = (a*yi +x[n-1])/(1.0+a) 
    for i in range(n-2,-1,-1):
        y[i] = yi = a*yi + b*y[i]

def mean(x):
    n = len(x) 
    return sum(x)/n

