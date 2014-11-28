# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 08:44:07 2014

@author: jarbona
"""
from lsimu import LSimu
from polymer import Polymer
from constrain import Box

def test_lammps_poly_simu():
    Simu = LSimu()
    Simu.create_ploymers(NP=10,N=100,box=Box([0,0,0],[10,10,10]))
    Simu.generate_xyz("test/mix.xyz",Mass="one")
    #Simu.generate_pdb("test/mix.pdb")
    #Simu.generate_script("test/basic.txt",run_size=10)
    #Simu.run
def test_lammps_poly_add_onep():
    Simu = LSimu()
    Simu.create_ploymers(NP=10,N=100,box=Box([0,0,0],[10,10,10]))
    Simu.add(Polymer(N=1))
    Simu.generate_xyz("test/mix.xyz",Mass="one")   
    Simu.generate_pdb("test/mix.pdb")
def test_lammps_from_hand():
    Simu = LSimu()
    box = Box([0,0,0],[10,10,10])
    P1 = Polymer(N=20,type_bead=1,ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[],gconstrain=[box])

    Simu.add(P1)
    Simu.add_box(box)
    Simu.generate_xyz("test/mix.xyz",Mass="one")
    Simu.generate_interactions("test/interactions")
    Simu.generate_pdb("test/mix.pdb")
    Simu.generate_script("test/basic.txt",run_length=1000,samplingrate=10,initconf="test/mix.xyz",
                         outtraj="test/out.dcd",outfile="test/out.xyz",interactions="test/interactions",particle="1")
    Simu.run(script="test/basic.txt")
    
def test_lammps_from_hand_mix():
    Simu = LSimu()
    box = Box([0,0,0],[10,10,10])
    P1 = Polymer(N=20,type_bead=1,liaison={"1-1":[1,1],"1-2":[1,2],"2-2":[1,3]},ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[],gconstrain=[box])
    P2 = Polymer(N=20,type_bead=2,liaison={"1-1":[1,1],"1-2":[1,2],"2-2":[1,3]},ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[],gconstrain=[box])

    Simu.add(P1)
    Simu.add(P2)
    Simu.add_box(box)
    Simu.generate_xyz("test/mix.xyz",Mass="one")
    Simu.generate_interactions("test/interactions")
    Simu.generate_pdb("test/mix.pdb")
    Simu.generate_script("test/basic.txt",run_length=1000,samplingrate=10,initconf="test/mix.xyz",
                         outtraj="test/out.dcd",outfile="test/out.xyz",interactions="test/interactions",particle="1 2")
    Simu.run(script="test/basic.txt")
    
def test_lammps_angle():
    Simu = LSimu()
    box = Box([0,0,0],[10,10,10])
    P1 = Polymer(N=20,type_bead=1,liaison={"1-1":[1,1]},angle_bond=True,angle_def={"1-1-1":[20,1]},ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[],gconstrain=[box])

    Simu.add(P1)
    Simu.add_box(box)
    Simu.generate_xyz("test/lp.xyz",Mass="one")
    Simu.generate_interactions("test/ainteractions")
    Simu.generate_pdb("test/lp.pdb")
    Simu.generate_script("test/lp.txt",run_length=10000,samplingrate=10,initconf="test/lp.xyz",
                         outtraj="test/lp.dcd",outfile="test/lp.xyz",interactions="test/ainteractions",particle="1 2")
    Simu.run(script="test/lp.txt")
    
def test_lammps_from_hand_mix2():
    Simu = LSimu()
    box = Box([0,0,0],[10,10,10])
    P1 = Polymer(N=20,type_bead=[1]*10+[2]*10,liaison={"1-1":[1,1],"1-2":[1,2],"2-2":[1,3]},ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[],gconstrain=[box])

    Simu.add(P1)
    Simu.add_box(box)
    Simu.generate_xyz("test/mix.xyz",Mass="one")
    Simu.generate_interactions("test/interactions")
    Simu.generate_pdb("test/mix.pdb")
    Simu.generate_script("test/basic.txt",run_length=1000,samplingrate=10,initconf="test/mix.xyz",
                         outtraj="test/out.dcd",outfile="test/out.xyz",interactions="test/interactions",particle="1 2")
    Simu.run(script="test/basic.txt")
    
def test_lammps_from_hand_mix_interactions():
    Simu = LSimu()
    box = Box([0,0,0],[10,10,10])
    P1 = Polymer(N=20,type_bead=[1]*10+[2]*10,liaison={"1-1":[1,1],"1-2":[1,2],"2-2":[1,3]},ptolerance=0,type_polymer="linear",start_id=0,lconstrain=[],gconstrain=[box])

    Simu.add(P1)
    Simu.add_box(box)
    Simu.generate_xyz("test/mix.xyz",Mass="one")
    Simu.generate_interactions("test/softinteractions")
    Simu.generate_pdb("test/mix.pdb")
    Simu.generate_script("test/basic.txt",run_length=1000,samplingrate=10,initconf="test/mix.xyz",outtraj="test/out.dcd",
                                     outfile="test/final.xyz",interactions="test/softinteractions",particle="1 2")
    Simu.run(script="test/basic.txt")
    
    box = Box([-5,-5,-5],[15,15,15])
    Simu.add_box(box)

    Simu.generate_xyz("test/mix2.xyz",Mass="one",from_lammps_xyz="test/final.xyz")
    Simu.generate_interactions("test/interactions",bond="fene")
    Simu.generate_script("test/basic_rest.txt",template_name="./template/basic.txt",
                         run_length=1000,samplingrate=10,initconf="test/mix2.xyz",outtraj="test/out2.dcd",outfile="test/out2.xyz",
                         interactions="test/interactions",particle="1 2")
    Simu.run(script="test/basic_rest.txt")
if __name__ == "__main__":
    test_lammps_from_hand_mix_interactions()