# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 17:35:41 2014

@author: jarbona
"""

from create import one_polymer
import numpy as np
from utils import norm
import unittest
from constrain import Box,Point
 
def test_create():
    
    coords,bonds,type_beads,ids = one_polymer(N=2,type_bead=0,liaison={"0-0":[1.0,0]},ptolerance=0,type_polymer="linear",start_id=0)
    assert (bonds[0] == [[0,0,0,1]])
    assert (round(norm(coords[0]-coords[1]),3) == 1.00)
    assert(type_beads == [0,0])
    
def test_create_mixture():
    coords,bonds,type_beads,ids = one_polymer(N=3,type_bead=[0,1,0],liaison={"0-0":[1.0,0],"0-1":[2.0,1]},ptolerance=0,type_polymer="linear",start_id=0)
    assert (bonds[0] == [[0,1,0,1],[1,1,1,2]])
    assert (round(norm(coords[0]-coords[1]),3) == 2.00)
    assert (round(norm(coords[1]-coords[2]),3)== 2.00)
def test_startid():   
    coords,bonds,type_beads,ids = one_polymer(N=2,type_bead=0,liaison={"0-0":[1.0,0]},ptolerance=0,type_polymer="linear",start_id=1)
    assert (bonds[0] == [[0,0,1,2]])
    
def test_one_particule():   
    coords,bonds,type_beads,ids = one_polymer(N=1,type_bead=0,liaison={"0-0":[1.0,0]},ptolerance=0,type_polymer="linear",start_id=1)
    assert (bonds[0] == [])

def test_startbond():   
    coords,bonds,type_beads,ids = one_polymer(N=2,type_bead=0,liaison={"0-0":[1.0,0]},ptolerance=0,type_polymer="linear",start_id=1,start_bond=1)
    assert (bonds[0] == [[1,0,1,2]])
    
def test_angle():   
    coords,bonds,type_beads,ids = one_polymer(N=2,type_bead=0,liaison={"0-0":[1.0,0]},ptolerance=0,type_polymer="linear",start_id=1,start_bond=1,angle_bond=True)
    assert (bonds[1] == [])
    
def test_create_mixture_angle():
    coords,bonds,type_beads,ids = one_polymer(N=3,type_bead=[0,1,0],
                                             liaison={"0-0":[1.0,0],"0-1":[2.0,1]},
                                              angle_bond=True,
                                              angle_def={"0-0-0":[1.,0],"0-0-1":[1.,1],"0-1-1":[1.,2],"1-1-1":[1.,3]},
                                        ptolerance=0,type_polymer="linear",start_id=0)
    assert (bonds[0] == [[0,1,0,1],[1,1,1,2]])
    assert (bonds[1] == [[0,1,0,1,2]])
def test_create_mixture_angle_start_angle():
    coords,bonds,type_beads,ids = one_polymer(N=3,type_bead=[0,1,0],
                                               liaison={"0-0":[1.0,0],"0-1":[2.0,1]},
                                              angle_bond=True,start_angle=10,
                                              angle_def={"0-0-0":[1.,0],"0-0-1":[1.,1],"0-1-1":[1.,2],"1-1-1":[1.,3]},
                                        ptolerance=0,type_polymer="linear",start_id=0)
    assert (bonds[0] == [[0,1,0,1],[1,1,1,2]])
    assert (bonds[1] == [[10,1,0,1,2]])

def test_bond_size():
    box = Box([-100,-100,-100],[400,400,400])


    coords,bonds,type_beads,ids = one_polymer(N=400,type_bead=1,liaison={"1-1":[2.0,1]},ptolerance=0,
                                       type_polymer="linear",start_id=0,lconstrain=[],
                                       gconstrain=[box],max_trial=300000,rc=0.5,virtual_lp=None,rigid_constrain=True)  
                                       
    bonds = np.mean(np.sqrt(np.sum((coords[1:]-coords[:-1])**2,axis=1)))
    assert(abs(bonds-2)< 0.1)
    
def test_bond_size_angle():
    box = Box([-100,-100,-100],[400,400,400])
    
    coords,bonds,type_beads,ids = one_polymer(N=400,type_bead=1,liaison={"1-1":[1.5,1]},ptolerance=0,
                                       type_polymer="linear",start_id=0,lconstrain=[],
                                       gconstrain=[box],max_trial=500000,rc=1,virtual_lp=4,rigid_constrain=False,flexible_lp=True)  
                                       
    bonds = np.mean(np.sqrt(np.sum((coords[1:]-coords[:-1])**2,axis=1)))
    assert(abs(bonds-1.5)< 0.1)                                       
                                       
def test_crash():
    PS = [0.2,0.2,0.2]
    PM = [10,10.2,0.2]
    PE = [0.2,199,19]
    start = Point(index=0,position=PS)  
    middle = Point(index=39,position=PM)
    
    end = Point(index=399,position=PE)
    #end = Point(index=199,position=[0.2,0.2,0.2])
    box = Box([-10,-10,-10],[400,400,400])
    coords_lp,bonds,type_beads,ids = one_polymer(N=400,type_bead=1,liaison={"1-1":[1.0,1]},ptolerance=0,
                                           type_polymer="linear",start_id=0,lconstrain=[start,middle,end],
                                           gconstrain=[box],max_trial=500000,rc=0.5,virtual_lp=3,rigid_constrain=True)

def test_create_mixture_angle_start_angle_start_id():
    coords,bonds,type_beads,ids = one_polymer(N=3,type_bead=[0,1,0],
                                               liaison={"0-0":[1.0,0],"0-1":[2.0,1]},
                                              angle_bond=True,start_angle=10,
                                             angle_def={"0-0-0":[1.,0],"0-0-1":[1.,1],"0-1-1":[1.,2],"1-1-1":[1.,3]},
                                        ptolerance=0,type_polymer="linear",start_id=20)
    assert (bonds[0] == [[0,1,20,21],[1,1,21,22]])
    assert (bonds[1] == [[10,1,20,21,22]])
    #Test starting id
    
def test_constraint():
    
    box = Box()
    assert (box.is_inside([0.5,0.5,0.5]))
    assert ( not box.is_inside([0.5,0.5,1.5]))
    assert ( not box.is_inside([0.5,1.5,0.5]))
    assert ( not box.is_inside([1.5,0.5,0.5]))
    
def test_constrain_box_on_polymer():
    box = Box([0,0,0],[10,10,10])
    coords,bonds,type_beads,ids = one_polymer(N=200,type_bead=1,ptolerance=0,type_polymer="linear",start_id=0,gconstrain=[box])
    for p in coords:
        assert(box.is_inside(p))
        
def test_constrain_on_polymer_start():    
    start = Point(index=0,position=[0.2,0.2,0.2])
    #end = Point(index=199,position=[0.2,0.2,0.2])
    box = Box([0,0,0],[10,10,10])
    coords,bonds,type_beads,ids = one_polymer(N=20,type_bead=1,ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[start],gconstrain=[box])   
    assert(norm(coords[0]-np.array([0.2,0.2,0.2])) < 0.01)
    
def test_constrain_on_polymer_quasi_start():    
    start = Point(index=1,position=[0.2,0.2,0.2])
    #end = Point(index=199,position=[0.2,0.2,0.2])
    box = Box([0,0,0],[20,20,20])
    coords,bonds,type_beads,ids = one_polymer(N=20,type_bead=1,ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[start],gconstrain=[box],max_trial=1000)   
    print norm(coords[0]-np.array([0.2,0.2,0.2]))
    assert(norm(coords[0]-np.array([0.2,0.2,0.2])) < 3)
    
def test_constrain_on_polymer_end():    
    end = Point(index=19,position=[0.2,0.2,0.2])
    #end = Point(index=199,position=[0.2,0.2,0.2])
    box = Box([0,0,0],[20,20,20])
    coords,bonds,type_beads,ids = one_polymer(N=20,type_bead=1,ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[end],gconstrain=[box],max_trial=1000)   
    print norm(coords[-2]-np.array([0.2,0.2,0.2]))
    assert(norm(coords[-2]-np.array([0.2,0.2,0.2])) < 3)
    
def test_constrain_on_polymer_start_end():
    start = Point(index=1,position=[0.2,0.2,0.2])    
    end = Point(index=199,position=[0.2,0.2,0.2])
    #end = Point(index=199,position=[0.2,0.2,0.2])
    box = Box([0,0,0],[20,20,20])
    coords,bonds,type_beads,ids = one_polymer(N=200,type_bead=1,ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[start,end],gconstrain=[box],max_trial=1000)   
    #print norm(coords[-1]-np.array([0.2,0.2,0.2]))
    print norm(coords[-2]-np.array([0.2,0.2,0.2]))
    print norm(coords[1]-np.array([0.2,0.2,0.2]))
    assert(norm(coords[-2]-np.array([0.2,0.2,0.2])) < 3)
    assert(norm(coords[1]-np.array([0.2,0.2,0.2])) < 3)

def test_constrain_on_polymer_start_middle_end():
    start = Point(index=1,position=[0.2,0.2,0.2])  
    middle = Point(index=100,position=[0.2,0.2,0.2])

    end = Point(index=199,position=[0.2,0.2,0.2])
    #end = Point(index=199,position=[0.2,0.2,0.2])
    box = Box([0,0,0],[20,20,20])
    coords,bonds,type_beads,ids = one_polymer(N=200,type_bead=1,ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[start,middle,end],gconstrain=[box],max_trial=1000)   
    #print norm(coords[-1]-np.array([0.2,0.2,0.2]))
    print norm(coords[-2]-np.array([0.2,0.2,0.2]))
    print norm(coords[1]-np.array([0.2,0.2,0.2]))
    assert(norm(coords[-2]-np.array([0.2,0.2,0.2])) < 3)
    assert(norm(coords[99]-np.array([0.2,0.2,0.2])) < 3)
    assert(norm(coords[101]-np.array([0.2,0.2,0.2])) < 3)
    assert(norm(coords[1]-np.array([0.2,0.2,0.2])) < 3)
    
class MyTestCase(unittest.TestCase):
    
    def testerror1(self):
        self.assertRaises(SomeCoolException, one_polymer(N=3,type_bead=[0,1],liaison={"0-0":[1.0,0]},ptolerance=0,type_polymer="linear",start_id=0))