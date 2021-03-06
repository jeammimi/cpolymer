Inferring the physical properties of yeast chromatin through Bayesian analysis of whole nucleus simulations
=======

 [![DOI](https://zenodo.org/badge/27133678.svg)](https://zenodo.org/badge/latestdoi/27133678)

This notebook explain how to create an initial configuration for the Yeast nucleus, how to set a persistence length
and how to run the simulation.
The simulations will generate a trajectory file and also dumps the contact to compute hic
matrices

http://nbviewer.ipython.org/github/jeammimi/cpolymer/blob/master/notebook/Yeast-Nucleus-persistence.ipynb



In the following we detail all the other possibilities of this library

cpolymer
========

Library to create polymer path with given functionnality, such as
some constrain on the position of any monomer.
It is also possible to create polymer with some rigidity

What can be achieved with it:

https://www.youtube.com/watch?v=iwep58wIsuU

To install
=======

The easiest way:
sudo pip install cpolymer

Example of uses
=======

Here I show how to use the module to creato a block copolymer system:

http://nbviewer.ipython.org/github/jeammimi/cpolymer/blob/master/notebook/cpolymer_demo.ipynb

Here How to create a polymer whose path has some specific constrain

http://nbviewer.ipython.org/github/jeammimi/cpolymer/blob/master/notebook/polymer_path.ipynb

Finally a little bit of fun with the creation of a polymeric eiffel tower

http://nbviewer.ipython.org/github/jeammimi/cpolymer/blob/master/notebook/cpolymer_eiffel_demo.ipynb

How to create a simple nucleus

http://nbviewer.ipython.org/github/jeammimi/cpolymer/blob/master/notebook/Simple-nucleus.ipynb

How to create a Yeast nucleus

http://nbviewer.ipython.org/github/jeammimi/cpolymer/blob/master/notebook/Yeast-Nucleus.ipynb

How to create a Yeast nucleus with a persistence length and dumping hic

http://nbviewer.ipython.org/github/jeammimi/cpolymer/blob/master/notebook/Yeast-Nucleus-persistence.ipynb

How to create a Yeast and make it replicate!

http://nbviewer.ipython.org/github/jeammimi/cpolymer/blob/master/notebook/Cell_division.ipynb


Here is the video from the last notebook!

https://www.youtube.com/watch?v=iwep58wIsuU


To visualise with vmd
=======
First install vmd then use the script  visu/rep.tcl :
for example:

vmd  -e ../visu/rep.tcl ./nucleus_yeastnoyau2.pdb dump_init.nucleus_yeast.comp.dcd final.xyz
