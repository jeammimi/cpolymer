# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 09:17:44 2014

@author: jarbona
"""
import numpy as np

class Constrain:
    def __init__(self):
        pass
    def is_inside(pt):
        return True
        
class Box(Constrain):
    def __init__(self,bl=[0,0,0],tr=[1,1,1]):
        Constrain.__init__(self)
        self.bl = bl 
        self.tr = tr
    def is_inside(self,pt):
        for x1,x2,x in zip(self.bl,self.tr,pt):
            if not(x1< x <x2):
                return False
        return True
    def generate(self):
        p = []
        for x1,x2 in zip(self.bl,self.tr):
            p.append((x2-x1)*np.random.rand()+x1)
        return p
class Point:
    def __init__(self,index=0,position=[0,0,0]):
        self.index = index
        self.position = position
        
        