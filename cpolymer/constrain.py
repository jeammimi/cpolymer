# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 09:17:44 2014

@author: jarbona
"""
import numpy as np
from .utils import norm, generateV
class Constrain:
    def __init__(self):
        pass
    def is_inside(pt):
        return True
        
class Box(Constrain):
    def __init__(self,bl=[0,0,0],tr=[1,1,1]):
        Constrain.__init__(self)
        self.bl = np.array(bl)
        self.tr = np.array(tr)
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
    def rescale(self,v=1):
        for i in range(3):
            self.bl[i] *=v
            self.tr[i] *= v
    @property
    def center(self):
        return self.tr/2 + self.bl/2
class Sphere(Constrain):
    def __init__(self,position=[0,0,0],radius=1):
        Constrain.__init__(self)
        self.position = position
        self.radius = np.array(radius)
    def is_inside(self,pt):
        if norm(self.position-np.array(pt)) <= self.radius:
            return True
        else:
            return False
    def generate(self):
        r = self.radius * np.random.random()
        return self.position + r*generateV()
    def rescale(self,v=1):
        for i in range(3):
            self.position[i] *=v
        self.radius *= v
class Point:
    def __init__(self,index=0,position=[0,0,0]):
        self.index = index
        self.position = position
    def __repr__(self):
        return "%i %s"%(self.index,str(self.position))
        
        