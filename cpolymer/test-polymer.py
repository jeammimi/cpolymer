# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 08:54:40 2014

@author: jarbona
"""

from polymer import Polymer
from constrain import Box
import numpy as np

def test_unwrap():
    
    l=10
    box = Box([0,0,0],[l,l,l])
    P1 = Polymer(N=20,type_bead=1,ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[],gconstrain=[box])
    
    for i in range(20):
        P1.coords[i][0] = i
        if P1.coords[i][0] > l:
            P1.coords[i][0] -= l
    print P1.coords
    P1.unwrap(ref=[0,0,0],box=box)
    print P1.coords
    assert np.sum(P1.coords[::,0] - np.array(range(20))) < 1e-9
    
def test_unwrap2():
    
    l=10
    box = Box([0,0,0],[l,l,l])
    P1 = Polymer(N=20,type_bead=1,ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[],gconstrain=[box])
    
    for i in range(20):
        P1.coords[i][0] = i
        if P1.coords[i][0] > l:
            P1.coords[i][0] -= l +0.5
            
    print P1.coords
    P1.unwrap(ref=[10,10,10],box=box)
    print P1.coords
    assert np.sum(P1.coords[::,0] - np.array(range(10,20)+[i+0.5 for i in range(20,30)])) < 1e-9