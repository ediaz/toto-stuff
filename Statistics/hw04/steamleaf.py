 
"""
steamleaf.py
 
Author:
   Ernesto P. Adorio, Ph.D.
  UPDEPP at Clark Field,
  Pampanga, the Philippines
"""
 
from   math import *
import scipy.stats as stat
 
 
def stemleafpairs(X, stempos= 0, leafwidth=1):
    """
    X - data array
    stempos   - position of last digit of stem,
                  from decimal point.
    leafwidth - number of digits in leaves.
 
    Return value:
     a list of stem-leaf pairs
    """
    eps = 10**-10
    stem10 = pow(10, stempos)
    leaf10 = 10**leafwidth
 
    output = []
    for x in X:
        y = x
        if stempos > 0:
           leaf, stem = modf (x * stem10)
        else:
           leaf, stem = modf(x/ stem10)
 
        leaf = abs(leaf * leaf10)
#        print x, int(stem), round(leaf) # decomment after testing!
        output.append((int(stem), int(leaf+eps) ))
    return output
 
 
def prettyprint(Pairs, stemwidth=4):
    """
    Given a list of Pairs (stem, leaf), prints it out for
    """
    Pairs.sort()
    minstem = Pairs[0][0]
    maxstem = Pairs[-1][0]
 
 
    laststem = minstem
    outstr = "%*d"  % (stemwidth, minstem*2) + "|"
    for (stem, leaf) in Pairs:
       if stem!= laststem:
          print outstr
          outstr = ""
          for i in range(laststem+1, stem):
             outstr = "%*d" % (stemwidth,  i*2) + "|"
             print outstr
          outstr =   "%*d" % (stemwidth, stem*2) + "|"
          laststem = stem
       outstr += str(leaf)
    print outstr
