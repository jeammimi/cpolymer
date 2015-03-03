# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 08:45:43 2014

@author: jarbona
"""
import numpy as np
import bisect


def init_proba(P,dx,lrange=[0,3.14]):
    """
    Generate a uni dimentionnal probability
    by creating a cumulative sum following the P law
    P must depand only on one variable x
    """
 
    x = np.arange(lrange[0],lrange[1],dx)

    Plist  = np.array([P(xx) * dx for xx in x])
    norm  = sum( Plist )
    Plist /= norm
    #print "Norm" , norm
    try:
        index = [Plist[0] ]
    except:
        print lrange[0],lrange[1],dx
    for w,xx in enumerate(x[1:]):
        index.append(index[-1] + Plist[w+1] )
    return x,index
    
def init_proba_log(P,dx,lrange=[0,3.14]):
    """
    Generate a uni dimentionnal probability
    by creating a cumulative sum following the P law
    P must depand only on one variable x
    """
 
    x = np.arange(lrange[0],lrange[1],dx)
    
    Plist  = np.array([P(xx)  for xx in x])
    #print Plist
    
    Plist -= max(Plist)
    Plist = np.exp(Plist)
    #print Plist
    
    norm  = sum( Plist )
    Plist /= norm
    #print "Norm" , norm
    try:
        index = [Plist[0] ]
    except:
        print lrange[0],lrange[1],dx
    for w,xx in enumerate(x[1:]):
        index.append(index[-1] + Plist[w+1] )
    return x,index
    
def generate_point_proba(x,index):
    try:
        u = np.random.rand()
        i = bisect.bisect(index,u)
        #print index[i-1],u,index[i] 
        return x[i]
    except:
        
        raise "Constrain to high"
    
