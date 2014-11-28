# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 09:19:38 2014

@author: jarbona
"""
import inspect
from create import one_polymer
class Polymer:
    def __init__(self,**kwargs):
        coords,bonds,type_beads,ids = one_polymer(**kwargs)
        
        #Getting liaison and angle information
        self.liaison = kwargs.get("liaison",None)
        self.angle_def = kwargs.get("angle_def",None)
        
         
        a = inspect.getargspec(one_polymer)
        for arg,value in zip(a.args[-len(a.defaults):],a.defaults):
            #print arg
            if self.liaison is None and arg == "liaison":
                 self.liaison =value
            if self.angle_def is None and arg == "angle_def":
                 self.angle_def = value
            
        self.coords = coords
        self.bond = bonds[0]
        self.angle = bonds[1]
        self.types_beads = type_beads
        self.ids = ids
    def change_coords(self,coords):
        if len(coords) == len(self.coords):
            self.coords = coords
        
    def get_types_beads(self):
        return set(self.types_beads)
        
    def get_xyz(self,start_id=1,start_bond=0,start_angle=0):#,start_dihedral=0):
        
        #Must set the id of the atom to the start_id
        id_add = start_id - self.ids[0]
        bond_add = start_bond - self.bond[0][0]
        angle_add = 0
        if self.angle != []:
            angle_add = start_angle - self.angle[0][0]
        #dihedral_add = start_bond - self.bond[0][0]
        Atom,Bond,Angle = [],[],[]
        
        for ids,pos,name in zip(self.ids,self.coords,self.types_beads):
            X,Y,Z = pos
            Atom.append("%10i%10i%10i%10.4f%10.4f%10.4f\n"%(ids+id_add, name,name,X,Y,Z))
        for bid,typebond,n1,n2 in self.bond:
            Bond.append("%10i%10i%10i%10i\n"%(bid + bond_add, typebond, n1 + id_add, n2 + id_add))
        for aid,typebond,n1,n2,n3 in self.angle:
            Angle.append("%10i%10i%10i%10i%10i\n"%(aid + angle_add, typebond,n1 + id_add,n2 + id_add,n3 + id_add))
        return Atom,Bond,Angle
        