# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 18:06:14 2014

@author: jarbona
"""
import random
import numpy as np
from math import sin,cos,sqrt

def generateV():
    theta = random.uniform(0,2*3.14)
    phi = random.uniform(0,3.14)
    return np.array([cos(theta) * sin(phi),sin(theta) * sin(phi),cos(phi)])
def generate_rotation_m(v,theta):
    x,y,z=v
    s = sin(theta)
    c = cos(theta)
    M = np.array([[x**2+(1-x**2) * c   , x * y *(1-c)-z * s , x *z * (1-c)+y * s ],
                      [x *y *(1-c)+ z * s , y**2+(1-y**2)* c , y * z * (1-c)-x * s ],
                       [x *z *(1-c)-y* s , y * z * (1-c)+x * s , z**2+(1-z**2)* c]])
    
    return M

def norm(v):
    n = 0
    for el in v:
        n += el *el
    return np.sqrt(n)
def generate_angle(v,angle=3.14*60/180.):
    """
    From an inital vector v 
    return a random vector rotated from angle
    """
    a = random.uniform(-0.3,0.3)
    b = random.uniform(-0.3,0.3)
    r = random.choice([0,1,2])
    if r == 2:
        vp = np.array([a,b,(-a*v[0]-b*v[1])/(v[2])])
    if r == 1:
        c=b
        vp = np.array([a,(-a*v[0]-c*v[2])/(v[1]),c])
    if r == 0:
        c=a
        vp = np.array([(-b*v[1]-c*v[2])/(v[0]),b,c])
    vp /= sqrt(sum(vp*vp))
    M = generate_rotation_m(vp,theta=angle)
    v = np.dot(M,v)
    return  v / np.sqrt(sum(v*v))