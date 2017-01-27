# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 17:57:56 2014

@author: jarbona
"""
import types
from .utils import generateV,norm
import numpy as np
from .proba import init_proba,init_proba_log,generate_point_proba

def one_polymer(N=2,type_bead=1,liaison={"1-1":[1.0,1]},type_polymer="linear",
                angle_bond=False,angle_def={"1-1-1":[1.0,1]},
                start_id=1,start_bond=0,start_angle=0,
                ptolerance=0,
                gconstrain=[],lconstrain=[],max_trial=300000,rc=0.5,virtual_lp=None,rigid_constrain=True,flexible_lp=True):
    """
    Function to generate a polymer chain
    return coordinate list , [bond list,angle list] , type of bead list, list of id of the monomers
    
    N is the total number of monomers
    type_bead can be set as an integer whose value us the type of monomer
              or with  a list containing the type of bead of each monomer
    liaison is a dictionnary which contain in key the name of the type of bead
            linked. The first value of the item is the size of the length and
            the second the type of bonds (Used later to defined ist properties with lammps)
            Must be key must be specified as 1-2 and not 2-1 for example
    angle_bond set to True to define bond between three monomer
    angle_def is the same as liaison but to define angle bond between three monomers
    start_id is used to shift the ids of the monomers
    start_bond is used to shift the ids of the bonds
    start_angle is used to shift the ids of the angle
    ptolerance (not implemented)
    gconstrain : list of global spatial constrained to be respected
    lconstrain : list of local contstrain to be respected
    max_trial : maximium number of trial to generate a point when some contrain are applied
    rc is used to enforce the constrain
    virtual_lp if different of none (Ex 2.0) the rigidity of the initial chain will be increased
    Be carefull the value here won t be proportionnal to a persistence length
    rigid_constrain if set to true and a local constrain is defined, the monomer will be exaclty at that position
    flexible_lp: if True and virtual_lp is not None: every max_trial/5 virtual_lp is divided by 1.5
    """
    
    assert(N > 0)
    
    ids = [i for i in range(N)]
    
    if type(type_bead) == int:
        type_beads = [type_bead for i in range(N)]
    else:
        type_beads = type_bead
         
    if type_polymer == "linear":
        #zero is set for the initialisation
        #the real type of bond is computed later
        bonds = [[start_bond + i,0,i, i+1] for i in range(N - 1)]
        angles = []
        if angle_bond:
            if N >=3:
                angles = [[start_angle + i,0,i,i+1,i+2] for i in range(N - 2)]
          
    if type_polymer == "circular":
        bonds = [[start_bond + i,0,i, i+1] for i in range(N )] 
        bonds[-1][-1] = 0
        angles = []
        if angle_bond:
            if N >=3:
                angles = [[start_angle + i,0,i,i+1,i+2] for i in range(N )]
                
                angles[-1][-2] = 1
                angles[-1][-1] = 0
                angles[-2][-1] = 0
                
    assert(len(type_beads) == N)
    
   
    bond_sizes = []
    for nb,bond in enumerate(bonds):
        idb,typebond,m1,m2 = bond
        tl = [type_beads[m1],type_beads[m2]]
        tl.sort()
        b1,b2 = tl
       
        bond_sizes.append(liaison["%i-%i"%(b1,b2)][0] )       
        bonds[nb][1] = liaison["%i-%i"%(b1,b2)][1]
    
    if angle_bond:
        angle_sizes = []
        for nb,angle in enumerate(angles):
            ida,typeangle,m1,m2,m3 = angle
            
            tl = [type_beads[m1],type_beads[m2],type_beads[m3]]
            tl.sort()
            b1,b2,b3 = tl
           
            #print b1,left,b2
            angle_sizes.append(angle_def["%i-%i-%i"%(b1,b2,b3)][0] )       
            angles[nb][1] = angle_def["%i-%i-%i"%(b1,b2,b3)][1]
    
    coords = []
    coords.append(generate_next(coords,gconstrain=gconstrain,lconstrain = lconstrain,max_trial=max_trial,bond_sizes=bond_sizes,rc=rc,
                                virtual_lp=virtual_lp,rigid_constrain=rigid_constrain,flexible_lp=flexible_lp))
    
    for bond,lbond in zip(bonds,bond_sizes):
        coords.append(generate_next(coords,gconstrain=gconstrain,lconstrain = lconstrain,max_trial=max_trial,bond_sizes=bond_sizes,rc=rc,
                                    virtual_lp=virtual_lp,rigid_constrain=rigid_constrain,flexible_lp=flexible_lp))
        
    if start_id != 0:
        bonds = [ [bond[0],bond[1],bond[2] + start_id,bond[3] + start_id] for bond in bonds ]
        angles = [ [angle[0],angle[1],angle[2] + start_id,angle[3] + start_id,angle[4] + start_id] for angle in angles ]

        ids = [idt + start_id for idt in ids]
    
    return np.array(coords) ,[bonds,angles],type_beads ,ids


def generate_from_local_constrain(coords,bond_sizes=[],lconstrain=[],max_trial=100,rc=0.1,virtual_lp=None,rigid_constrain=True):
    pos = []
    
    index_point = len(coords)
    
    # To disregard the constrain that where before the actual index we have to 
    # find the first constrain that is usefull
    # The constrain should be organized
    #print [0 if c.index < index_point else 1 for c in lconstrain ]
    is_there_constrain = [0 if c.index < index_point else 1 for c in lconstrain ]
    #print virtual_lp
    if (1 not in is_there_constrain) and (virtual_lp is None):
        #No more constrain
        start = np.array(coords[-1])
        pos = start + bond_sizes[index_point-1] * generateV()
        return np.array(pos),[]
   
    else:
        if 1 in is_there_constrain:
            start_constrain = is_there_constrain.index(1)
        else:
            start_constrain= None
        
    if start_constrain is not None and lconstrain[start_constrain].index == index_point and rigid_constrain:
        # if a constrain match the index we return the position
        return np.array(lconstrain[start_constrain].position),[]
    
    redo = []
    for xi in range(3):
        
        
        extent = []
        width = []
        if start_constrain != None:
            for c in lconstrain[start_constrain:]:
                extent.append(c.position[xi])
                width.append(np.sum(bond_sizes[index_point:c.index]))
                if not rigid_constrain and width[-1] == 0:
                    width[-1]=rc
                elif width[-1] == 0:
                    width[-1] = 1
            #print width, "h"
        if coords != []:
            extent.append(coords[-1][xi])
            width.append(bond_sizes[index_point-1])
        #print width,bond_sizes[:10]
            
        if len(coords) >= 2 and virtual_lp:
            #width[-1] is the bond length
            extent.append(coords[-1][xi] + width[-1]*(coords[-1][xi]-coords[-2][xi]) / norm(coords[-1]-coords[-2]))
            width.append(1./virtual_lp)

        
        width_quad = np.array(width) #** 2  * rc
        cm = 1/np.sum(1/width_quad) * np.sum(np.array(extent)/np.array(width_quad))
        
        gauss = lambda x,c,w : np.exp(-(x-c)**2/(2*w)) #Not square as we use width_quad
        def Proba(x):
           
            prod = []
            for pos,w in zip(extent,width_quad):
                
                prod.append(gauss(x,pos,w=w))
            return np.prod(prod)
            
        lngauss = lambda x,c,w : -(x-c)**2/(2*w)
        def lnProba(x):
           
            prod = []
            for pos,w in zip(extent,width_quad):
                
                prod.append(lngauss(x,pos,w=w))
            return np.sum(prod)
        w = np.min(width)
        w = max(1.,w)
        #print w,cm,extent
        #print width
        x,index = init_proba_log(lnProba,dx=w/5.,lrange=[cm-8*w,cm+8*w])
        #x,index = init_proba(Proba,dx=w/5.,lrange=[cm-4*w,cm+4*w])
        #print bond_sizes[index_point-1:]
        p = generate_point_proba(x,index)
        
        #print width ,width_quad
        #if coords != []:
        #    print extent , cm ,coords[-1][xi]
        #exit()
        
#        if xi == 0:
#            #print 
#            print index_point , start_constrain 
#            print extent
#            #print width
#            if len(coords) >= 2 :
#                print coords[-1],2*coords[-1]-coords[-2], "CM",cm , p
#            else:
#                print "CM",cm , p
            #print p
        pos.append(p)
        redo.append([x,index])
    return np.array(pos),redo
    
def generate_next(coords,gconstrain=[],lconstrain=[],max_trial=100,bond_size=1,bond_sizes=[],rc=0.1,virtual_lp=None,rigid_constrain=True,flexible_lp=True):
    
    N = 0
    redo = []  # The probability densitiy is generated once and then several trial are extracted from it
    virtual_lp0 = virtual_lp  #virtual_lp0 will be decrease if flexible_lp in True
    while N < max_trial:
        
        if coords == [] and lconstrain == []:
            if gconstrain == []:
                pos =  np.array([0,0,0])
            else:
                pos = np.array(gconstrain[0].generate())
              
        else:
            if flexible_lp and virtual_lp != None and N != 0 and N % int(max_trial/5) == 0:
                virtual_lp0 /= 2.
                redo = []
                print(("Deacreasing virtual_lp to ", virtual_lp0))
                
            if redo == []:
                
                pos,redo = generate_from_local_constrain(coords,lconstrain=lconstrain,bond_sizes=bond_sizes,rc=rc,virtual_lp=virtual_lp0,rigid_constrain=rigid_constrain)
            else:
                pos  = np.array( [ generate_point_proba(x,index) for x,index in redo ])
    
        # check global constrain
             
        out=False
        for gc in gconstrain:
            if not gc.is_inside(pos):
                out = True
        if out:
            N += 1
                    
            if N == max_trial-1:
                print(coords)
                print("constrain not satisfied")
                raise 
            continue 
        break
    return pos
  

    