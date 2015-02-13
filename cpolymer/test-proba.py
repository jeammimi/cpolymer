# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 09:00:06 2014

@author: jarbona
"""
import numpy as np
from proba import init_proba,generate_point_proba,init_proba_log

def test_proba():
    p = lambda x,c,w : np.exp(-(x-c)**2/(2*w**2))
    
    x,index = init_proba(P=lambda x : p(x,1,2),dx=0.1,lrange=[-10,10])
    
    pts = np.array([generate_point_proba(x,index) for i in range(100000)])
    
    print pts.mean()
    print pts.std()
    assert(abs(pts.mean()-1) < 0.1)
    assert(abs(pts.std()-2) < 0.1)
    
def test_proba_log():
    
    p = lambda x,c,w : -(x-c)**2/(2*w**2)
    
    x,index = init_proba_log(P=lambda x : p(x,1,2),dx=0.1,lrange=[-10,10])
    #print index
    #print init_proba(P=lambda x : np.exp(p(x,1,2)),dx=0.11,lrange=[-10,10])[1]
    
    pts = np.array([generate_point_proba(x,index) for i in range(100000)])
    
    print pts.mean()
    print pts.std()
    assert(abs(pts.mean()-1) < 0.1)
    assert(abs(pts.std()-2) < 0.1)