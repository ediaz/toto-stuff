import numpy as np




def mix_2gauss(n,m1,m2,s1,s2,alpha):
    if alpha >1 : alpha=1

    U = uniform(0,1,n)

    I=(U<alpha)

    g1=np.random.normal(m1,s1,n)
    g2=np.random.normal(m2,s2,n)

    mix = I*g1+(1-I)*g2 
   
 
    return mix 


def gauss(m,s,n):
    g1=np.random.normal(m,s,n)
    return g1

def uniform(a,b,n):
    u=np.random.uniform(a,b,n)

    return u

def lognormal(m,s,n):
    ln = np.random.lognormal(m,s,n)
    return ln


def histogram_plot(x,n,fig=None):
    import matplotlib.pylab as plt
    from   matplotlib.pyplot import  savefig
    
    count,bins,ignored = plt.hist(x,n, normed=True)
    plt.xlim([-3,100])
    plt.ylim([0,.7])
    if fig==None:
        plt.show()
    else:
        savefig(fig,dpi=200,bbox_inches=0)
 
