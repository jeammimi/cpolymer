# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 09:19:38 2014

@author: jarbona
"""
import inspect
from create import one_polymer
import numpy as np
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
        self.extrabond = []
    def change_coords(self,coords):
        if len(coords) == len(self.coords):
            self.coords = coords
    @property        
    def n_beads(self):
        return len(self.types_beads)

        
    def shift(self,add_id=0,add_bond=0,add_angle=0):
        self.ids = [iid + add_id for iid in self.ids]
        self.bond = [ [b[0]+add_bond,b[1],b[2] + add_id,b[3] + add_id ] for b in self.bond]
        self.angle = [ [b[0]+add_angle,b[1],b[2] + add_id,b[3] + add_id,b[4] + add_id ] for b in self.angle]
        
    def get_types_beads(self):
        return set(self.types_beads)
    
    def get_atoms_bonds_angle(self,start_id=1,start_bond=0,start_angle=0):
        id_add = start_id - self.ids[0]
        bond_add = 0
        if self.bond != []:
            bond_add = start_bond - self.bond[0][0]
        angle_add = 0
        if self.angle != []:
            angle_add = start_angle - self.angle[0][0]
        #dihedral_add = start_bond - self.bond[0][0]
        Atoms,Bonds,Angles = [],[],[]
        
        for ids,pos,name in zip(self.ids,self.coords,self.types_beads):
            X,Y,Z = pos
            Atoms.append([ids+id_add, name,X,Y,Z])
        for bid,typebond,n1,n2 in self.bond:
            Bonds.append([bid + bond_add, typebond, n1 + id_add, n2 + id_add])
        for aid,typebond,n1,n2,n3 in self.angle:
            Angles.append([aid + angle_add, typebond,n1 + id_add,n2 + id_add,n3 + id_add])
        return Atoms,Bonds,Angles
    
    def get_xyz_format_lammps(self,start_id=1,start_bond=0,start_angle=0):#,start_dihedral=0):
        
        #Must set the id of the atom to the start_id
        id_add = start_id - self.ids[0]
        bond_add = 0
        if self.bond != []:
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
    
    def bond_size(self):
        return np.sqrt(np.sum((self.coords[1:]-self.coords[:-1])**2,axis=1))
        
    def unwrap(self,ref,box):
        for i in range(self.n_beads):
            for xi in range(3):
                lxi = (box.tr[xi]-box.bl[xi])/2
                #print self.coords[i][xi]-ref[xi] , lxi
                if self.coords[i][xi]-ref[xi] > lxi:
                    n = round((self.coords[i][xi]-ref[xi])/(2*lxi))
                    #print "N",n
                    self.coords[i][xi] -= 2*n*lxi
                elif self.coords[i][xi]-ref[xi] < -lxi:
                    n = abs(round((self.coords[i][xi]-ref[xi])/(2*lxi)))
                    #print "N",n
                    self.coords[i][xi] += 2*n*lxi
            ref = self.coords[i]