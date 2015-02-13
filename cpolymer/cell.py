# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 14:27:43 2014

@author: jarbona
"""
from halley.constrain import Spherical,Nowhere
from halley.vectors import V
import numpy as np
from constrain import Sphere,Point
from polymer import Polymer

def get_ref_points(snucleus,smicrotubule,l_len_p=[]):
    center = Spherical(V.O,radius=snucleus)
    mt = Spherical(V(-snucleus,0,0),radius=smicrotubule)
    depart = Spherical(V(0,0,0),radius=snucleus-smicrotubule+0.5)
    I = depart*mt
    res_start = []
    for chromosome in l_len_p:
        if chromosome is None:
            res_start.append(None)
            continue
        m = I.get_random()
        res=[m]
        for x in chromosome:
            #print m,x,smicrotubule,snucleus
            if x > 1.5*snucleus:
                x = 1.5*snucleus
            sp = Spherical(m,radius=x)
            radius_x = center * sp
            #print isinstance(radius_x,Nowhere)
            if isinstance(radius_x,Nowhere):
                print "Warning no intersection for this size of polymer"
                if x < smicrotubule:
                    print "increasing size"
                    sp = Spherical(m,radius=smicrotubule)
                if x > snucleus:
                    print "decreasing size"
                    sp = Spherical(m,radius=1.5*snucleus)
                radius_x = center * sp
                
            start_x = radius_x.get_random()
            res.append(start_x)
        #Reorganise res
        res[0],res[1] = res[1],res[0]
        res_start.append(res)
    return res_start
    
def cell(chrlen,chrcen,Ribopos,lp,Radius,mt,liaison,angle_def,angle_bond):
    """
    chrlen and chrcen are in monomere unit
    Ribopos also
    lp also
    """
    chl = []
    chn = []
    for l,cen in zip(chrlen,chrcen):
        if cen == None:
            chl.append(None)
            chn.append([None,None,l])
        else:
            fact = 16
            
            chl.append([np.sqrt(2*fact*lp**2*cen),np.sqrt(2*fact*lp**2*(l-cen))])
            chn.append([cen,(l-cen),l])
    chl = np.array(chl)
  
    Refs = get_ref_points(snucleus=Radius,smicrotubule=mt,l_len_p=chl)
    #print liaison

    nucleus = Sphere(position=[0,0,0],radius=Radius)

    list_polymer = []
    for chromosome,(rep,length) in enumerate(zip(Refs,chn)):
        #print chromosome,length,chl[chromosome]
        #,rep
        #print Radius,mt
        
        #print start,middle,end
        type_bead = [1 for n in range(length[2])]
        if rep != None:
            #dangling chromosome    
            type_bead[0] = type_bead[-1] = 2
            type_bead[length[0]] = 4
            
        for ch,insert,rlength in Ribopos:
            if chromosome == ch-1:

                type_bead[insert:insert+rlength] = [3 for _i in range(rlength)]
                cm_nuc = insert + int(rlength)/2
                
        if rep != None:
            start = Point(index=0,position=rep[0]._v)
            middle = Point(index=length[0],position=rep[1]._v)
            end = Point(index=length[2],position=rep[2]._v)
            lconstrain=[start,middle,end]
        else:
            lconstrain = []
            
        for ch,insert,rlength in Ribopos:
            if chromosome == ch-1:

                nucleole = Point(index=cm_nuc,position=(0.66*Radius,0,0))
                lconstrain=[start,middle,nucleole,end]
        
        #print liaison
        list_polymer.append(Polymer(N=length[2],type_bead=type_bead,liaison=liaison,
                                                  angle_bond=angle_bond,
                                                  angle_def=angle_def,
                                                  ptolerance=0,type_polymer="linear",
                                                  lconstrain=lconstrain,gconstrain=[nucleus],
                                                  max_trial=300000,rc=1.,virtual_lp=None))
                                                  
    
    return list_polymer
if __name__ == "__main__":
    chrlen=[230218,
            813184,
            316620,
            1531933,
            576874,
            270161,
            1090940,
            562643,
            439888,
            745751,
            666816,
            1078177,
            924431,
            784333,
            1091291,
            948066 ]
            
    chrcen=[
            151465,
            238207,
            114385,
            449711,
            151987,
            148510,
            496920,
            105586,
            355629,
            436307,
            440129,
            150828,
            268031,
            628758,
            326584,
            555957]
    
                
    Ribopos=[(12,450000,750000)]
    
    C = 110
    w = 15
    lp = 15
    liaison={"0-0":1.0,"0-1":1.0,"0-2":1.0,"0-3":3,
            "1-1":1.0,"1-2":1,"1-3":3,
            "2-2":1.0,"2-3":3,
            "3-3":6},