# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 08:44:07 2014

@author: jarbona
"""
from cpolymer.lsimu import LSimu
from cpolymer.polymer import Polymer
from cpolymer.constrain import Box


def test_generate_interaction():
    Simu = LSimu()
    Simu.create_polymers(NP=10, N=100, box=Box([0, 0, 0], [10, 10, 10]))
    Simu.add_bond(typeb="harmonic", idbond=1, K=80, R0=1)
    Simu.add_pair(typep="lj/cut", idpair1=1, idpair2=1, epsilon=1, sigma=1, cutoff1=1.12)

    Simu.generate_interactions("test/interaction")


def test_lammps_poly_simu():
    Simu = LSimu()
    Simu.create_polymers(NP=10, N=100, box=Box([0, 0, 0], [10, 10, 10]))
    Simu.generate_xyz("test/mix.xyz", Mass="one")
    # Simu.generate_pdb("test/mix.pdb")
    # Simu.generate_script("test/basic.txt",run_size=10)
    # Simu.run


def test_lammps_poly_add_onep():
    Simu = LSimu()
    Simu.create_polymers(NP=10, N=100, box=Box([0, 0, 0], [10, 10, 10]))
    Simu.add(Polymer(N=1))
    Simu.generate_xyz("test/mix.xyz", Mass="one")
    Simu.generate_pdb("test/mix.pdb")


def test_lammps_from_hand():
    Simu = LSimu()
    box = Box([0, 0, 0], [10, 10, 10])
    P1 = Polymer(N=20, type_bead=1, ptolerance=0, type_polymer="linear",
                 start_id=0, lconstrain=[], gconstrain=[box])

    Simu.add(P1)
    Simu.add_box(box)
    Simu.add_bond(typeb="harmonic", idbond=1, K=80, R0=1)
    Simu.add_pair(typep="lj/cut", idpair1=1, idpair2=1, epsilon=1, sigma=1, cutoff1=1.12)

    Simu.generate_xyz("test/mixh.xyz", Mass="one")
    Simu.generate_interactions("test/interactionsh")
    Simu.generate_pdb("test/mixh.pdb")
    Simu.generate_script("test/basich.txt", run_length=1000, samplingrate=10, initconf="test/mixh.xyz",
                         outtraj="test/out.dcd", outfile="test/out.xyz", interactions="test/interactionsh", particle="1")
    Simu.run(script="test/basich.txt")


def test_lammps_from_hand_mix():
    Simu = LSimu()
    box = Box([0, 0, 0], [10, 10, 10])
    liaison = {"1-1": [1, 1], "1-2": [1, 2], "2-2": [1, 3]}
    P1 = Polymer(N=20, type_bead=1, liaison=liaison, ptolerance=0,
                 type_polymer="linear", start_id=0, lconstrain=[], gconstrain=[box])
    P2 = Polymer(N=20, type_bead=2, liaison=liaison, ptolerance=0,
                 type_polymer="linear", start_id=0, lconstrain=[], gconstrain=[box])

    Simu.add(P1)
    Simu.add(P2)
    for idl, value in list(liaison.items()):
        R0, idbond = value
        idpair1, idpair2 = list(map(int, idl.split("-")))
        Simu.add_bond(typeb="harmonic", idbond=idbond, K=80, R0=R0)
        Simu.add_pair(typep="lj/cut", idpair1=idpair1, idpair2=idpair2,
                      epsilon=1, sigma=R0, cutoff1=1.15)
    Simu.add_box(box)
    Simu.generate_xyz("test/mix.xyz", Mass="one")
    Simu.generate_interactions("test/interactionsm")
    Simu.generate_pdb("test/mix.pdb")
    Simu.generate_script("test/basic.txt", run_length=1000, samplingrate=10, initconf="test/mix.xyz",
                         outtraj="test/out.dcd", outfile="test/out.xyz", interactions="test/interactionsm", particle="1 2")
    Simu.run(script="test/basic.txt")


def test_lammps_angle():
    Simu = LSimu()
    box = Box([0, 0, 0], [10, 10, 10])
    liaison = {"1-1": [1, 1]}
    angle_def = {"1-1-1": [20, 1]}
    P1 = Polymer(N=20, type_bead=1, liaison=liaison, angle_bond=True, angle_def=angle_def,
                 ptolerance=0, type_polymer="linear", start_id=0, lconstrain=[], gconstrain=[box])

    Simu.add(P1)
    Simu.add_box(box)

    for idl, value in list(liaison.items()):
        R0, idbond = value
        idpair1, idpair2 = list(map(int, idl.split("-")))
        Simu.add_bond(typeb="harmonic", idbond=idbond, K=80, R0=R0)
        Simu.add_pair(typep="lj/cut", idpair1=idpair1, idpair2=idpair2,
                      epsilon=1, sigma=R0, cutoff1=1.15)
    for idl, value in list(angle_def.items()):
        K, idangle = value
        idpair1, idpair2, idpair3 = list(map(int, idl.split("-")))
        Simu.add_angle(typea="harmonic", idangle=idangle, K=K, theta=180)

    Simu.generate_xyz("test/lp.xyz", Mass="one")
    Simu.generate_interactions("test/ainteractions")
    Simu.generate_pdb("test/lp.pdb")
    Simu.generate_script("test/lp.txt", run_length=10000, samplingrate=10, initconf="test/lp.xyz",
                         outtraj="test/lp.dcd", outfile="test/lp.xyz", interactions="test/ainteractions", particle="1 2")
    Simu.run(script="test/lp.txt")


def test_lammps_from_hand_mix2():
    Simu = LSimu()
    box = Box([0, 0, 0], [10, 10, 10])
    liaison = {"1-1": [1, 1], "1-2": [1, 2], "2-2": [1, 3]}
    P1 = Polymer(N=20, type_bead=[1] * 5 + [2] * 5 + [1] * 5 + [2] * 5, liaison=liaison,
                 ptolerance=0, type_polymer="linear", start_id=0, lconstrain=[], gconstrain=[box])

    Simu.add(P1)

    Simu.add_box(box)
    for idl, value in list(liaison.items()):
        R0, idbond = value
        idpair1, idpair2 = list(map(int, idl.split("-")))
        Simu.add_bond(typeb="harmonic", idbond=idbond, K=80, R0=R0)
        if idpair1 == idpair2:
            cutoff1 = 3
        else:
            cutoff1 = 1.5
        Simu.add_pair(typep="lj/cut", idpair1=idpair1, idpair2=idpair2,
                      epsilon=1, sigma=R0, cutoff1=cutoff1)
    Simu.generate_xyz("test/mix.xyz", Mass="one")
    Simu.generate_interactions("test/interactions")
    Simu.generate_pdb("test/mix.pdb")
    Simu.generate_script("test/basic.txt", run_length=1000, samplingrate=10, initconf="test/mix.xyz",
                         outtraj="test/out.dcd", outfile="test/out.xyz", interactions="test/interactions", particle="1 2")
    Simu.run(script="test/basic.txt")


def test_lammps_from_hand_mix_interactions():
    Simu = LSimu()
    box = Box([0, 0, 0], [10, 10, 10])
    liaison = {"1-1": [1, 1], "1-2": [1, 2], "2-2": [1, 3]}
    P1 = Polymer(N=20, type_bead=[1] * 10 + [2] * 10, liaison=liaison, ptolerance=0,
                 type_polymer="linear", start_id=0, lconstrain=[], gconstrain=[box])

    Simu.add(P1)
    Simu.add_box(box)
    for idl, value in list(liaison.items()):
        R0, idbond = value
        idpair1, idpair2 = list(map(int, idl.split("-")))
        Simu.add_bond(typeb="harmonic", idbond=idbond, K=350, R0=R0)
        Simu.add_pair(typep="lj/cut", idpair1=idpair1, idpair2=idpair2,
                      epsilon=1, sigma=R0, cutoff1=1.15)

    Simu.generate_xyz("test/mix.xyz", Mass="one")
    Simu.generate_interactions("test/softinteractions")
    Simu.generate_pdb("test/mix.pdb")
    Simu.generate_script("test/basic.txt", run_length=10000, samplingrate=100, initconf="test/mix.xyz", outtraj="test/out.dcd",
                         outfile="test/final.xyz", interactions="test/softinteractions", particle="1 2")
    Simu.run(script="test/basic.txt")

    box = Box([0, 0, 0], [10, 10, 10])
    Simu.add_box(box)
    Simu.clean_interactions()
    for idl, value in list(liaison.items()):
        R0, idbond = value
        idpair1, idpair2 = list(map(int, idl.split("-")))
        Simu.add_bond(typeb="fene", idbond=idbond, K=80, R0=R0, epsilon=1.2, sigma=1)
        Simu.add_pair(typep="lj/cut", idpair1=idpair1, idpair2=idpair2,
                      epsilon=1, sigma=R0, cutoff1=1.15)

    Simu.generate_xyz("test/mix2.xyz", Mass="one", from_lammps_xyz="test/final.xyz")
    Simu.generate_interactions("test/interactions")
    Simu.generate_script("test/basic_rest.txt", template_name="./template/basic.txt",
                         run_length=1000, samplingrate=10, initconf="test/mix2.xyz", outtraj="test/out2.dcd", outfile="test/out2.xyz",
                         interactions="test/interactions", particle="1 2")
    Simu.run(script="test/basic_rest.txt")
if __name__ == "__main__":
    test_lammps_from_hand_mix_interactions()
