# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 09:42:32 2014

@author: jarbona
"""
import subprocess
import string
import copy
import numpy as np
from polymer import Polymer
from bond import Bond
from pair import Pair
from angle import Angle
import os

from sortn import sort_nicely


class LSimu:
    def __init__(self,cmd="lammps"):
        self.cmd=cmd
        self.molecules = []
        self.script = None
        self.xyz_name = None
        self.pdb_name = None
        self.box = None
        self.Bond = []
        self.Angle = []
        self.Pair = []
        self.extra_bond =[]
        self.level = 0
        #Generate list to define interaction for lammps simulation
        self.clean_interactions()
        
    def coarse(self,avoid_type=[],root="./cg",nlevel=2,periodic=True,lastcmd="lammps",**args):
        #Coarse_two_by two:
        self.Coarse = LSimu(cmd=self.cmd)
        self.Coarse.level = self.level + 1
         
        self.Coarse.box = copy.deepcopy(self.box)
        self.Coarse.box.rescale(1/2**0.5)
        Nbead = 0 
        for nmol,molecule in enumerate(self.molecules):
            #Coarse by type:
            nt = [molecule.types_beads[0]]
            translate=[0]
            for ix,c in enumerate(molecule.types_beads[1:],1):
                same_type = ( c  == nt[-1] )
                
                correspond_to_only_one_bead = not(len(translate) >= 2 and translate[-1] == translate[-2])
                #print same_type,correspond_to_only_one_bead,translate,nt
                if same_type and correspond_to_only_one_bead and c not in avoid_type:
                    translate.append(len(nt) -1)
                else:
                    nt.append(c)
                    translate.append(len(nt) - 1)
                    
            #Compute new CM
            if periodic:
                molecule.unwrap(ref=self.box.center,box=self.box)
            CM = np.zeros((len(nt),3))
            N = np.zeros(len(nt))
            for b in translate:
                CM[b] += molecule.coords[b]
                N[b] += 1
            angle_bond = False
            if molecule.angle != []:
                angle_bond = True
            Nbead += len(nt)
            P = Polymer(N=len(nt),type_bead=nt,gconstrain=[self.Coarse.box],liaison=molecule.liaison,angle_def=molecule.angle_def,angle_bond=angle_bond)
            #P.change_coords(CM/N[::,np.newaxis]/2**0.5)
            self.molecules[nmol].translate = [] + translate
            self.Coarse.add(P)
        print("At level {0}, {1} segments".format(self.Coarse.level,Nbead))
        #Take care of the extra bonds
        newb = []
        for mol1,mol2,typeb in self.extra_bond:
            
            nbt = [[mol1[0],self.molecules[mol1[0]-1].translate[mol1[1]] - 1],
                                  [mol2[0],self.molecules[mol2[0]-1].translate[mol2[1]] - 1],
                                     typeb]
            #When coarse graining extra bond can attach same molecules
            if nbt not in newb:
                newb.append(nbt)
                self.Coarse.add_extra_bond(nbt[0],nbt[1],nbt[2])
                                     
        #Try to take care of the interactions
        for b in self.iBond:
            self.Coarse.iBond.append(copy.deepcopy(b))
        for p in self.iPair:
            self.Coarse.iPair.append(copy.deepcopy(p))
            if p.typep == "lj/cut":
                self.Coarse.iPair[-1].args["sigma"] *= 2**(0.33-0.5)
                if "cutoff1" in self.Coarse.iPair[-1].args:
                    self.Coarse.iPair[-1].args["cutoff1"] *= 2**(0.33-0.5)

            if p.typep == "soft":
                self.Coarse.iPair[-1].args["A"] *= 2**(0.33-0.5)
        for a in self.iAngle:
            self.Coarse.iAngle.append(copy.deepcopy(a))
            self.Coarse.iAngle[-1].args["K"] /= 2**0.5
        
       
        
        self.Coarse.Mass = self.Mass
        self.Coarse.info_bond = self.info_bond
        self.Coarse.info_pair = self.info_pair
        self.Coarse.natom = self.natom
        
        rep = "{0}_{1}".format(root,self.Coarse.level)
        if not os.path.exists(rep):
            os.mkdir(rep)
            
            
        if self.Coarse.level != nlevel:
            self.Coarse.coarse(avoid_type=avoid_type,root=root,nlevel=nlevel,**args)
        #else:
        #    self.Coarse.cmd=lastcmd
            
        #Generate_xyz
        initconf = rep+"/in.xyz"
       
        
        
        self.Coarse.generate_xyz(initconf,Mass=self.Mass)
        #Generate_interactions
        pdb = rep +"/in.pdb"
        self.Coarse.generate_pdb(pdb)

        interactions = rep+"/interactions"
        self.Coarse.generate_interactions(interactions,info_bond = self.info_bond,info_pair=self.info_pair)
        #Generate_script
        outfile = rep+"/out.xyz"
        outtraj= rep + "/out.dcd"
        script = rep+"/cg.txt"
        run_length = args.get("run_length",1000)
        samplingrate = args.get("samplingrate",100)
        template_name = args.get("template_name","./template/basic.txt")
        self.Coarse.generate_script(script_name=script,template_name=template_name,
                                                        run_length=run_length,
                                                        samplingrate=samplingrate,
                                                        initconf=initconf,outtraj=outtraj,
                                                        outfile=outfile,
                                                        interactions=interactions,
                                                        particle=" ".join(map(str,range(1,max(self.natom)+1))),
                                                        **args)
        #Run
        self.Coarse.run(script=script)
        
        self.reassign_from_coarse(outfile,periodic=periodic)
        
    def assign_from_lammps(self,name,unwrap=False):
        coords = self.read_from_lammps(name)
        
        start = 0
        for molc,molecule in enumerate(self.molecules):
            self.molecules[molc].coords = coords[start:start+len(molecule.coords)]
            if unwrap:
                #print self.box.tr
                self.molecules[molc].unwrap(ref=self.box.center,box=self.box)
            start += len(molecule.coords)
                
    def reassign_from_coarse(self,outfile,periodic):
        
        self.Coarse.assign_from_lammps(outfile,unwrap=periodic)
        
         
            
        for nmol,(molecule,cmolecule) in enumerate(zip(self.molecules,self.Coarse.molecules)):
            #Get reverse mapping
            
            cmolecule.rtranslate = [[] for i in range(len(cmolecule.types_beads))]
            
            for i in range(len(molecule.types_beads)):
                cmolecule.rtranslate[molecule.translate[i]].append(i)
            #print cmolecule.rtranslate
            #Then assign position
            for molc,mapping in enumerate(cmolecule.rtranslate):
                if len(mapping) == 1:
                    self.molecules[nmol].coords[mapping[0]] = cmolecule.coords[molc]
                elif len(mapping) == 2:
                    if len(cmolecule.types_beads) == 1:
                        self.molecules[nmol].coords[mapping[0]] = cmolecule[molc].coords +np.random.rand(3)
                        self.molecules[nmol].coords[mapping[1]] = cmolecule[molc].coords +np.random.rand(3)
                    else:
                        if molc == 0:
                            t = cmolecule.coords[molc+1]-cmolecule.coords[molc]
                            l =  np.sqrt(np.sum(t**2))
                            if l == 0:
                                l = 1
                            t /= l
                            t0=t
                            t1 = t
                        elif molc == len(cmolecule.types_beads) - 1:
                            t = cmolecule.coords[molc]-cmolecule.coords[molc-1]
                            l =  np.sqrt(np.sum(t**2))
                            if l == 0:
                                l = 1
                            t /= l
                            t0 = t
                            t1=t
                        else:
                            t1 = cmolecule.coords[molc+1]-cmolecule.coords[molc]
                            l =  np.sqrt(np.sum(t1**2))
                            if l == 0:
                                l = 1
                            t1 /= l
                            t0 = cmolecule.coords[molc]-cmolecule.coords[molc-1]
                            l =  np.sqrt(np.sum(t0**2))
                            if l == 0:
                                l = 1
                            t0 /= l
                            
                        #print molc, l , cmolecule.coords[molc] - 0.3*l*t0 , cmolecule.coords[molc] + 0.3*l*t1
                        self.molecules[nmol].coords[mapping[0]] = cmolecule.coords[molc] - 0.25*l*t0
                        self.molecules[nmol].coords[mapping[1]] = cmolecule.coords[molc] + 0.25*l*t1
                            
                else:
                    print("something went wrong")
                    
            self.molecules[nmol].coords *= 2**0.5
            
            #Check it is inside the box
            if not periodic:
                for p in self.molecules[nmol].coords:
                    if not self.box.is_inside(p):
                        print "Mol not in the box"
                        raise
            else:
                #Should rewrap?
                pass
    def add(self,molecule):
        """
        add a molecule and take care of shifting the ids bonds
        """
        
        if self.molecules == []:
            start_id = 1
            start_bond = 0
            start_angle = 0
        else:
            start_id = self.molecules[-1].ids[-1] + 1
            
            #look for the firts previous molecule with a bond
            for mol in self.molecules[::-1]:
                if mol.bond != []:
                    start_bond = mol.bond[-1][0] + 1
                    
            if  self.molecules[-1].angle != []:
                start_angle = self.molecules[-1].angle[-1][0]
                
        add_id = start_id - molecule.ids[0]
        add_bond = 0
        if molecule.bond != []:
            add_bond = start_bond - molecule.bond[0][0]
        add_angle = 0
        if molecule.angle != []:
            add_angle = start_angle - molecule.angle[0][0]
        molecule.shift(add_id,add_bond,add_angle)
        self.molecules.append(molecule)
        
        
        
    def add_bond(self,**bond):
        """
        Adding bond interaction information:
        Sim = LSimu()
        Sim.add_bond(typeb="harmonic",idbond=1,K=60,R0=1.))
        Sim.add_bond(typeb="fene",idbond=2,K=80,R0=R0,epsilon=1.2,sigma=1)
        """
        self.iBond.append(Bond(**bond))
        
    def add_pair(self,**pair):
        self.iPair.append(Pair(**pair))
        
    def add_angle(self,**angle):
        self.iAngle.append(Angle(**angle))
    def add_box(self,box):
        self.box = box
        
    def create_polymers(self,NP,**kwargs):
        """ 
        create NP polymers
        All the option to cleate a polymer are defined in the function
        one_polymer in the create module
        if a box is specified it will be added in the global constrained of the polymer
        """
        if  kwargs.has_key("box"):
            self.box = kwargs["box"]
            kwargs.pop("box")
            if kwargs.has_key("gconstrain"):
                
                kwargs["gconstrain"].append(self.box)
            else:
                kwargs["gconstrain"] = [self.box]
        else:
            print("must specify a box") 
            raise
        
        self.molecules = []
        for p in range(NP):
            self.add(Polymer(**kwargs))
        
    def read_from_lammps(self,from_lammps_xyz):
        coords = []
        with open(from_lammps_xyz,"r") as f:
            for line in f:
                
                if line.startswith("Atoms"):
                    for line in f:
                        if line != "\n":
                            coords.append(map(float,line.split()[:6]))
                        else:
                            if coords != []:
                                break
                if coords != []:
                    break
          
            coords.sort()
            coords = np.array(coords)
            coords = coords[:,3:6]
        return coords
    def generate_xyz(self,xyz_name,Mass=None,from_lammps_xyz=None):
        """
        Generate xyz file to be read by lammps as an input file
        It is generated according to the molecules given to the LSimu object
        
        Mass can be set to one if all particle
        have a mass of one
        
        from_lammps_xyz can be specified. It should be the path of a file
        generated by lammps (by the command write_data)
        This file should contain the same number of particule than the system
        defined for the simulation. It is read and the initial position of
        the generated xyz file will be the same as the from_lammps_xyz file

        """
        self.xyz_name =xyz_name        
        
        start_id = 1
        start_bond = 0
        start_angle = 0
        self.Atom, self.Bond, self.Angle = [],[],[]
        coords = []
        if from_lammps_xyz is not None:
            coords = self.read_from_lammps(from_lammps_xyz)
           
        self.liaison = None
        self.angle_def = None
        self.natom = set()
        
        for molecule in self.molecules:
            if coords != []:
                molecule.change_coords(coords[start_id - 1:start_id - 1 + len(molecule.ids)])
            #still assigning a start_id but it should already be correct
            Atom,Bond,Angle = molecule.get_xyz(start_id=start_id,start_bond=start_bond,start_angle=start_angle)
            start_id += len(Atom) 
            start_bond += len(Bond)
            start_angle += len(Angle)
            self.Atom.extend(Atom)
            self.Bond.extend(Bond)
            self.Angle.extend(Angle)
            self.natom = set(self.natom).union(molecule.get_types_beads())
            
            if self.liaison is None:
                self.liaison = molecule.liaison
            else:
                if self.liaison != molecule.liaison:
                    #try to merge
                    pass
            
            if self.angle_def is None:
                self.angle_def = molecule.angle_def
            else:
                if self.angle_def != molecule.angle_def:
                    #try to merge
                    pass
                
        #Extra bonds:
        for mol1,mol2,typeb in self.extra_bond:
            absolute1 = self.molecules[mol1[0]-1].ids[mol1[1]-1]
            absolute2 = self.molecules[mol2[0]-1].ids[mol2[1]-1]
            bid = len(self.Bond)
            self.Bond.append("%10i%10i%10i%10i\n"%(bid, typeb, absolute1 , absolute2 ))
               
        
        
        self.ntype_bond = len(self.liaison.keys())
        
        self.ntype_angle = 0
        if self.angle_def is not None and self.Angle != []:
            self.ntype_angle = len(self.angle_def.keys())
        
        f = open(self.xyz_name,'w')
        f.write("\n")
        f.write("%8i     atoms\n"%len(self.Atom))
        f.write("%8i     bonds\n"%len(self.Bond)  )
        f.write("%8i     angles\n"%len(self.Angle))
        f.write("%8i     dihedrals\n"%(0))
        
        f.write("""
                 %i     atom types
                 %i     bond types
                 %i     angle types
                 0     dihedral types"""%(max(self.natom),self.ntype_bond,self.ntype_angle))
        if self.box is None:
            raise ("Must add a box: "
                    "ex: LSimu.add_box([0,0,0],[10,10,10])")
            
            
        f.write("""         
            %8.4f   %8.4f xlo xhi         
            %8.4f   %8.4f ylo yhi
            %8.4f   %8.4f zlo zhi
        
        Masses\n    
        """%(self.box.bl[0],self.box.tr[0],self.box.bl[1],self.box.tr[1],self.box.bl[2],self.box.tr[1]))
 
        if Mass == "one":
            Mass = { "%i"%i:1 for i in range(1,max(self.natom)+1)}
        elif Mass == None:
            print "Must set mass to one to give all atoms mass one or specify a mass for each atoms"
            raise
        self.Mass = Mass
        f.write("\n".join([ "         %i          %.1f"%(i,Mass["%i"%(i)]) for i in range(1,max(self.natom)+1)  ]))
        f.write("\nAtoms\n\n%s\n"%("".join(self.Atom)))
        f.write("Bonds\n\n%s\n"%("".join(self.Bond)))
        if self.Angle != []:
            f.write("Angles\n\n%s\n"%("".join(self.Angle)))
        f.close()
        
    def add_extra_bond(self,mol1,mol2,typeb):
        """
        The references start from one to follow lammps convention
        mol1 = [molid,atomid]
        mol2 = [molid,atomid]
        """
        self.extra_bond.append([mol1,mol2,typeb])
    
    def generate_pdb(self,pdb_name,shift=0,traduction={"1":"bead","2":"telo","3":"ribo","4":"cent","5":"spbb","6":"rcut","7":"scut"}):
        """
        Generate a pdb file according to the molecules given to the LSimu object
        traduction must contain as key the type of bead and as item the pdb name
        (should be 4 letters name)
        
        (Very ugly code)
        """
        #Initial stuff for pdb:
        towrite = []
        
        atomid = 0
       
        for n,molecule in enumerate(self.molecules):

            for rn,((x,y,z),type_bead) in enumerate(zip(molecule.coords,molecule.types_beads)):
                atomid+= 1
                
                name="bead"
                residue_name="bea"
                chain_id=string.ascii_letters[(n + shift) % len(string.ascii_letters)]
                residue_number=rn

                name = traduction["%i" % type_bead]
                
        
                towrite.append(self.one_pdb_atom(x,y,z,atomid,name,residue_name,residue_number,chain_id))
                
        with open(pdb_name,"w") as f:
            f.write("".join(towrite))
    def one_pdb_atom(self,x,y,z,atomid,name,residue_name,residue_number,chain_id):
          
        alternate=" "
        code_insertion_residue=" "
        occupancy=1
        temp=0.0
        seg_id="    "
        element_symbole=" "
        charge=" "
        atom="ATOM  "
        return "%6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"%(atom,atomid, name, alternate,residue_name,
            			 chain_id,residue_number,code_insertion_residue,x,y,z,occupancy,temp,
            				seg_id,element_symbole,charge)
        
        return
        
    def clean_interactions(self):
        """
        reset laamps interactions
        """
        self.iBond = []
        self.iPair = []
        self.iAngle = []
        
        
        
    def generate_interactions(self,interaction_name,info_bond=[],info_pair=[],write=True,sort=True):
        
        self.info_bond = info_bond
        self.info_pair = info_pair
        typeb = set([bond.typeb for bond in self.iBond])
        hybrid =""
        if len(typeb) > 1:
            hybrid="hybrid"
        self.Bond_interaction = ["bond_style {0} {1}\n".format(hybrid," ".join(typeb))]
        
        if info_bond != []:
            self.Bond_interaction.extend(info_bond)
            
        Bond_interaction =[]
        for i in range(len(self.iBond)):
            if len(typeb) > 1:
                self.iBond[i].hybrid=True
            Bond_interaction.append(str(self.iBond[i]))
        if sort:
            sort_nicely(Bond_interaction)
        self.Bond_interaction.extend(Bond_interaction)
            
        
                
        typep = set([pair.typep for pair in self.iPair])
        hybrid =""
        if len(typep) > 1:
            hybrid="hybrid"
        self.Pair_interaction = ["pair_style {0} {1}\n".format(hybrid," ".join(typep))]
        
        
        if len(typep) == 1 and list(typep)[0] == "lj/cut":
            #print "nthnthnt"
            #raise
            maxcut = max([pair.args.get("cutoff1",0) for pair in self.iPair ])
            self.Pair_interaction = ["pair_style lj/cut {0:.2f}".format(maxcut)]
            self.Pair_interaction.append("pair_modify shift yes")
            
        if len(typep) == 1 and list(typep)[0] == "soft":
            maxcut = max([pair.args.get("cutoff",0) for pair in self.iPair ])
            self.Pair_interaction = ["pair_style soft {0:.2f}".format(maxcut)]
           
            
            
        self.Pair_interaction.extend(info_pair)
        
        Pair_interaction = []
        for i in range(len(self.iPair)):
            if len(typep) > 1:
                self.iPair[i].hybrid=True
            Pair_interaction.append(str(self.iPair[i]))
        if sort:
            sort_nicely(Pair_interaction)
        self.Pair_interaction.extend(Pair_interaction)
        
        
        typea = set([angle.typea for angle in self.iAngle])
        
        hybrid =""
        if len(typea) > 1:
            hybrid="hybrid"
        self.Angle_interaction = ["angle_style {0} {1}\n".format(hybrid," ".join(typea))]

        Angle_interaction = []
        for i in range(len(self.iAngle)):
            if len(typep) > 1:
                self.iAngle[i].hybrid=True
            Angle_interaction.append(str(self.iAngle[i]))
        if sort:
            sort_nicely(Angle_interaction)
        self.Angle_interaction.extend(Angle_interaction)
        
        if write:
            self.write_interactions(interaction_name)
            
    def write_interactions(self,interaction_name):
        with open(interaction_name,"w") as g:
            g.write("\n".join(self.Bond_interaction)+"\n\n")
            g.write("\n".join(self.Pair_interaction)+"\n\n")
            if self.Angle != []:
                g.write("\n".join(self.Angle_interaction))
        
        
        
    def generate_interactions_o(self,interaction_name,bond="harmonic",SPB=False,radius=10,cutoff=1.15,reducedv_factor=1,khun=1.):

        if bond == "harmonic" :
            temp1 = "bond_style      harmonic\n"
            temp2 = "bond_coeff %i %.2f %.2f"
            ene = 350
        if bond == "fene":
            if SPB:
                temp1 = "bond_style hybrid harmonic fene\n"
                temp2 = "bond_coeff %i fene %.2f %.2f %.2f %.2f"
            else:
                temp1 = "bond_style fene\n"
            
                temp2 = "bond_coeff %i %.2f %.2f %.2f %.2f"
            ene = 30 * 1
        Bond = [temp1]
        if bond == "fene":
            Bond.append("special_bonds fene\n")
        ene_ratio=1#.35
        #Bond.append("special_bonds 0.0 1.0 1.0\n")
        if SPB and False:
            #print radius
            Pair = ["pair_style  hybrid lj/cut 3.0 gauss/cut  %.2f \n"%(2*radius) + "pair_modify shift yes\n"]
        else:
            Pair = ["pair_style   lj/cut 1.4 \n" + "pair_modify shift yes\n"]
        Angle = ["angle_style cosine/delta\n"] 
        Angle = ["angle_style harmonic\n"]
        keyl = range(1,len(self.natom)+1)
        for t1 in keyl:
            for t2 in keyl:
                if t2 >= t1:
                    if not self.liaison.has_key("%s-%s"%(t1,t2)):
                        print "Warning liaison between {0} and {1} not defined".format(t1,t2)
                    dist,tybe_b = self.liaison["%s-%s"%(t1,t2)]
                    if  cutoff is not None:                   
                        cut = dist*cutoff
                    else:
                        cut = dist*pow(2.,1/6.) 
                    odist = copy.deepcopy(dist)
                    #dist=1
                    if bond == "fene":
                        Bond.append(temp2 % (tybe_b,ene_ratio * ene/(dist*dist),1.5*dist,ene_ratio,dist) +"\n")
                    else:
                        Bond.append(temp2 % (tybe_b, ene,dist) +"\n")
                    dist = odist
                    precise =""
                    if SPB and False:
                        precise = "lj/cut" 
                    reduced = 1
                    if t1 == t2 and t1 == 3:
                        reduced = 1
                    else:
                        reduced = reducedv_factor 
                        
            
                    Pair.append("""pair_coeff	 %s %s %s %.1f %.2f  %.2f\n"""%(t1,t2,precise,ene_ratio,dist*reduced,cut*reduced))
                    if self.angle_def is not None and self.Angle != []:
                        for t3 in keyl:
                            if t3 >= t2:
                                dist,tybe_b = self.angle_def["%s-%s-%s"%(t1,t2,t3)]
                                k=khun/2. * dist # dist = 1 if no ribo involved else 0
                                Angle.append("angle_coeff %i %.3f 180.0\n"%(tybe_b,k))
        if SPB:
            if bond == "fene":
                Bond.append("bond_coeff %i harmonic 0 0 \n"%(self.liaison["spb"][1]))
                spbond = "bond_coeff %i harmonic %.1f %.1f\n"%(self.liaison["spb"][1],10,microtubule/realsigma)
            else:
                Bond.append("bond_coeff %i  0 0 \n"%(self.liaison["spb"][1]))
                spbond = "bond_coeff %i %.1f %.1f\n"%(self.liaison["spb"][1],10,microtubule/realsigma)
            n_i = len(diameter)/2
            for t1 in range(len(diameter) ):
                Pair.append("""pair_coeff	 %i %i %s 0. %.2f  %.2f\n"""%(t1+1,num_particle["spb"],precise,dist,cut))
              
            Pair.append("""pair_coeff	 %i %i %s 0. %.2f  %.2f\n"""%(num_particle["spb"],num_particle["spb"],precise,dist,cut))
        
        g = open(interaction_name,"w")
        g.write("".join(Bond)+"\n")
        g.write("".join(Pair)+"\n")
        if self.Angle != []:
            g.write("".join(Angle))
    
    def generate_script(self,script_name,template_name="./template/basic.txt",**kwargs):
        
        
        #First look in the template directory
        if os.path.exists(os.path.join(os.path.dirname(__file__), template_name)):
            template_name = os.path.abspath(os.path.join(os.path.dirname(__file__), template_name))
            print "Reading template" , template_name
       
  
       
        with open(template_name,"r") as f:
            template = string.Template("".join(f.readlines()))
            template = template.safe_substitute(kwargs)
            
        with open(script_name,"w") as f:
            f.write(template)
        
    def run(self,script):
        self.script = script
        if self.script != None:
            print "{0} < {1}".format(self.cmd,self.script)
            output = subprocess.check_output("{0} < {1}".format(self.cmd,self.script), shell=True)
            return output
        else:
            print("No script found\n")
            raise