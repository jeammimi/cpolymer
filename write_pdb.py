# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 09:39:06 2013

@author: jarbona
"""

import numpy,random,math,copy,string, time
from math import cos ,sin
from translate_cdf_cdc import write_single_pdb
random.seed(time.time())
import sys


def generateV():
    theta = random.uniform(0,2*3.14)
    phi = random.uniform(0,3.14)
    return numpy.array([cos(theta) * sin(phi),sin(theta) * sin(phi),cos(phi)])
def generate_rotation_m(v,theta):
    x,y,z=v
    s = sin(theta)
    c = cos(theta)
    M = numpy.array([[x**2+(1-x**2) * c   , x * y *(1-c)-z * s , x *z * (1-c)+y * s ],
                      [x *y *(1-c)+ z * s , y**2+(1-y**2)* c , y * z * (1-c)-x * s ],
                       [x *z *(1-c)-y* s , y * z * (1-c)+x * s , z**2+(1-z**2)* c]])
    
    return M

def norm(v):
    n = 0
    for el in v:
        n += el *el
    return numpy.sqrt(n)
def generate_angle(v,angle=3.14*60/180.):
    """
    From an inital vector v 
    return a random vector rotated from angle
    """
    a = random.uniform(-0.3,0.3)
    b = random.uniform(-0.3,0.3)
    r = random.choice([0,1,2])
    if r == 2:
        vp = numpy.array([a,b,(-a*v[0]-b*v[1])/(v[2])])
    if r == 1:
        c=b
        vp = numpy.array([a,(-a*v[0]-c*v[2])/(v[1]),c])
    if r == 0:
        c=a
        vp = numpy.array([(-b*v[1]-c*v[2])/(v[0]),b,c])
    vp /= math.sqrt(sum(vp*vp))
    M = generate_rotation_m(vp,theta=angle)
    v = numpy.dot(M,v)
    return  v / numpy.sqrt(sum(v*v))

def bound(v,l):
    for d in v:
        if abs(d) > l:
            return False
    return True
    
def duplicate(chromosomes,centrolist,Ribopos,diameter,num_particle,Box=None,cell="yeast",angle_bond=False,SPB=False,cutoff=1.4,khun=None,Extra=[],cutlp=10,coords=[]):
    keyl = diameter.keys()
    n_i = len(diameter)
    for k in keyl:
        diameter[str(int(k) + n_i)] = diameter[k]
    chromosomes += chromosomes
    centrolist += centrolist
    
    #Duplicate ribopos
    rp = copy.deepcopy(Ribopos)
    for ch,start,length in rp:
        Ribopos.append([ch+n_i,start,length])
    #Duplicate extra
    if Extra != []:
        e = copy.deepcopy(Extra)
        for ch,start,length,typeb in e:           
            Extra.append([ch+n_i,start,length,typeb+n_i])
    coord,maxi = get_coord(chromosomes,centrolist,Ribopos,diameter,num_particle,Box=Box,cell=cell,angle_bond=angle_bond,
                           SPB=SPB,cutoff=cutoff,khun=khun,Extra=Extra,cutlp=cutlp,DUPLICATE=True)
    return coord,maxi
    
    
    
def get_coord(chromosomes,centrolist,Ribopos,diameter,Mass,num_particle,Box=None,cell="yeast",angle_bond=False,SPB=False,cutoff=1.4,khun=None,Extra=[],cutlp=10,DUPLICATE=False,microtubule=0.3,REP=""
            ,dists=None,boxs=None,startconf=None,reducedv_factor=1):
                
    """
    if startconf generetae inital configuration from startconfile
    """


 
        

    rigidcut=6
    
    #create a list of special element to put ribo telo centro...
    special = [[[0,num_particle["telo"]],[c,num_particle["centro"]]] if c != None else [] for l,c in zip(chromosomes,centrolist)]
    if cell == "plasmo":
        from var import Var
        for name,chromo,pos1,pos2 in Var:
            where = int(pos1/5.)
            if int(pos1/5.) >=chromosomes[chromo - 1]:
                print chromo,name, int(pos1/5.),chromosomes[chromo - 1]
                print "Var outside chromo"
                where = chromosomes[chromo - 1] -2
            special[chromo - 1].append([where+0,num_particle["telo"]])
        
   
    for ch,start,length in Ribopos:
        chromosomes[ch-1] += length
        for p in range(start,start+length):
            special[ch - 1].append([p,num_particle["ribo"]])
    

    if Extra != []:
        for ch,start,length,typeb in Extra:           
            chromosomes[ch-1] += length
            for p in range(start,start+length):
                special[ch - 1].append([p,typeb])
     
    for nc,(l,c) in enumerate(zip(chromosomes,centrolist)):
        if c != None:
            special[nc].append([l-1,num_particle["telo"]])
    Atomtype = len(diameter)


    lbond = {}
    langle = {}
    keyl = diameter.keys()
    
    Typebond = 0
    Anglebond = 0
    
    for key1 in keyl:
        for key2 in keyl:
            if int(key1) <= int(key2):
                Typebond += 1
                lbond[key1+key2] = [diameter[key1]/2.+diameter[key2]/2.,Typebond]
                lbond[key2+key1] = lbond[key1+key2]
                
                if angle_bond:
                    for key3 in keyl:
                        if int(key2) <= int(key3):
                            Anglebond += 1
                            if int(key1) == num_particle["ribo"] or int(key2) == num_particle["ribo"] or int(key3)== num_particle["ribo"]:
                                langle[key1+key2+key3] = [0,Anglebond]
                            elif int(key1) == rigidcut or int(key2) == rigidcut or int(key3)==rigidcut:
                                langle[key1+key2+key3] =  [ cutlp * (diameter[key1]/3.+diameter[key2]/3.+ diameter[key3]/3.),Anglebond]
                            else:
                                langle[key1+key2+key3] =  [diameter[key1]/3.+diameter[key2]/3.+ diameter[key3]/3.,Anglebond]
               
    if SPB:
        Typebond += 1
        lbond["spb"]=[True,Typebond]
        
    Sum=1
    Atom=[]
    Bond=[]
    Angle=[]
    bond=1
    angle=1
    maxi=0
    mini=100
    chlist = []
    namelist = []
    centrolist = []
    
    boule = 50.
    if dists is not None:
        boule = 0
   
    
    start = generateV()
    pos_start = start *random.uniform(boule/8,boule)
    vstart = start
     #1u
    angle_bead = 80

    coord = []
    for n,l in enumerate(chromosomes):
        chlist.append([[],[],[],[]])
        namelist.append([])
        start = generateV()
        pos_start = start *random.uniform(boule/8.,boule)
        vstart = start
        print "Start", start
        nameprevious=1
        namepreviousprevious=1
        for i in range(l):
            number=i + Sum
            name=1
    
            for sp in special[n]:
                if i == sp[0]:
                    name=sp[1]
                    
            if name == 4: #centro
                centrolist.append(number)
        #    print s+s+s+s+s+s+s+s
            if i != 0:
                #Bond section
                typebond = lbond["%i%i"%(name,nameprevious)][1]
                dist = 1.2 * lbond["%i%i"%(name,nameprevious)][0]
                Bond.append("%10i%10i%10i%10i\n"%(bond, typebond,number-1,number))
                bond += 1
                if i!= 1 and angle_bond:
                    la = [name,nameprevious,namepreviousprevious]
                    la.sort()
                    typebond = langle["%i%i%i"%(la[0],la[1],la[2])][1]
                    Angle.append("%10i%10i%10i%10i%10i\n"%(angle, typebond,number-2,number-1,number))
                    angle += 1
                    
                #small circular
                if i == l-1 and l == 3:
                    Bond.append("%10i%10i%10i%10i\n"%(bond, typebond,number-2,number))
                    bond += 1
                
            
            if i == 0:
                dist = 0
            #print dist, name,nameprevious
            if dists is not None:
                dist = dists
            vstart = generate_angle(vstart,3.14*angle_bead/180.)
            if Box:               
                while not bound(pos_start +vstart * dist,Box):
                    
                    vstart = generate_angle(vstart,3.14*angle_bead/180.)
            pos_start += vstart * dist
            X,Y,Z = copy.deepcopy(pos_start)
            if startconf:
                X,Y,Z = startconf[n][i]
            chlist[-1][0].append(X)
            chlist[-1][1].append(Y)
            chlist[-1][2].append(Z)            
            namelist[-1].append(name)
          
            maxi=max(maxi,abs(X),abs(Y),abs(Z))

            coord.append([X,Y,Z])
            Atom.append("%10i%10i%10i%10.4f%10.4f%10.4f\n"%(number, name,name,X,Y,Z))
            namepreviousprevious = copy.deepcopy(nameprevious)
            nameprevious = copy.deepcopy(name)
        Sum += l
    
    if SPB:
        Sum
    #    Typebond += 1
        Atomtype += 1
        spbp = -radius
        Atom.append("%10i%10i%10i%10.4f%10.4f%10.4f\n"%(Sum ,num_particle["spb"],num_particle["spb"],spbp,0,0))
        print  spbp
    #    maxi  = -spbp
    #    mini = spbp
        print mini,maxi
    #    for i,atom in enumerate(centrol):
    #        print sum( (centroC[i]-numpy.array([-10,0,0])) *  (centroC[i]-numpy.array([-10,0,0])) )
    #        Bond.append("%10i%10i%10i%10i\n"%(bond, Typebond,Sum,atom))
    #        bond += 1
            
        Sum += 1
    
    extraatom=""
    if SPB:
        alternate=" "
        code_insertion_residue=" "
        occupancy=1
        temp=0.0
        atom="ATOM  "
        seg_id="    "
        element_symbole=" "
        charge=" "
        number =Sum-1
        name="spbb"
        residue_name="bea"
        chain_id=string.ascii_letters[len(chlist) + 1]
        residue_number=188 #au pif
        X,Y,Z = -10.025,0.843,0.744 
        extraatom = "%6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"%(atom,number, name, alternate,residue_name,
        			 chain_id,residue_number,code_insertion_residue,X,Y,Z,occupancy,temp,
        				seg_id,element_symbole,charge)

        #extraatom = "ATOM   %i spbb bea %s 188     -10.025   0.843   0.744  1.00  0.00\n"%(Sum-1,string.ascii_letters[len(chlist) + 1])
        
        
    write_single_pdb(REP+"/%snoyau2.pdb"%cell,chlist,namelist,extraatom=extraatom)
    maxi *=2
    f = open(REP+"%sconf2.txt"%cell,'w')
    f.write("\n")
    f.write("%8i     atoms\n"%len(Atom))
    
    if SPB:
        for cent in centrolist:
            Bond.append("%10i%10i%10i%10i\n"%(bond,lbond["spb"][1],cent,Sum-1))
            bond += 1
        
    if boxs is not None:
        maxi = radius * 1.1
    else:
        maxi+10
    f.write("%8i     bonds\n"%len(Bond)  )
    f.write("%8i     angles\n"%len(Angle))
    f.write("%8i     dihedrals\n"%(0))
    f.write("""
             %i     atom types
             %i     bond types
             %i     angle types
             0     dihedral types"""%(Atomtype,Typebond,Anglebond))
    f.write("""         
        %8.4f   %8.4f xlo xhi         
        %8.4f   %8.4f ylo yhi
        %8.4f   %8.4f zlo zhi
    
    Masses\n    
"""%(-maxi,maxi,-maxi,maxi,-maxi,maxi))
    print diameter
    key = diameter.keys()
    key.sort()
    #f.write("\n".join([ "         %s          %.1f"%(i,diameter[i]**3)  for i in key  ]))
    f.write("\n".join([ "         %i          %.1f"%(i,Mass["%i"%(i)]) for i in range(1,Atomtype+1)  ]))
    f.write("\nAtoms\n\n%s\n"%("".join(Atom)))
    f.write("Bonds\n\n%s\n"%("".join(Bond)))
    if angle_bond:
        f.write("Angles\n\n%s\n"%("".join(Angle)))
    f.close()
    #ATP file 
    
    
    bond = generate_lammps(diameter,lbond,num_particle,langle=langle,DUPLICATE=DUPLICATE,SPB=SPB,REP=REP,reducedv_factor=reducedv_factor)
    softbond = generate_lammps(diameter,lbond,num_particle,langle=langle,DUPLICATE=DUPLICATE,SPB=SPB,REP=REP,reducedv_factor=reducedv_factor,bond="harmonic",name="softinteractions")

    return coord ,maxi /2. + 10,bond,softbond
    

    
def generate_lammps(diameter,lbond,num_particle,langle=[],name="interactions",DUPLICATE=False,SPB=False,REP="",reducedv_factor=1,bond = "fene"):
    
    spbond = ""
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
    keyl = diameter.keys()
    for t1 in keyl:
        for t2 in keyl:
            if t2 >= t1:
                dist,tybe_b = lbond["%s%s"%(t1,t2)]
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
                if angle_bond:
                    for t3 in keyl:
                        if t3 >= t2:
                            dist,tybe_b = langle["%s%s%s"%(t1,t2,t3)]
                            k=khun/2. * dist # dist = 1 if no ribo involved else 0
                            Angle.append("angle_coeff %i %.3f 180.0\n"%(tybe_b,k))


    if SPB:
        if bond == "fene":
            Bond.append("bond_coeff %i harmonic 0 0 \n"%(lbond["spb"][1]))
            spbond = "bond_coeff %i harmonic %.1f %.1f\n"%(lbond["spb"][1],10,microtubule/realsigma)
        else:
            Bond.append("bond_coeff %i  0 0 \n"%(lbond["spb"][1]))
            spbond = "bond_coeff %i %.1f %.1f\n"%(lbond["spb"][1],10,microtubule/realsigma)
        n_i = len(diameter)/2
        for t1 in range(len(diameter) ):
            Pair.append("""pair_coeff	 %i %i %s 0. %.2f  %.2f\n"""%(t1+1,num_particle["spb"],precise,dist,cut))
          
        Pair.append("""pair_coeff	 %i %i %s 0. %.2f  %.2f\n"""%(num_particle["spb"],num_particle["spb"],precise,dist,cut))
        
    if SPB and False:
        #Bond.append("""bond_coeff %i harmonic %.2f %.2f\n"""%(Typebond,350,3.0))
        n_i = len(diameter)/2
        for t1 in range(len(diameter) -1):
            Pair.append("""pair_coeff	 %i %i %s 0. %.2f  %.2f\n"""%(t1+1,num_particle["spb"],precise,dist,cut))
            if DUPLICATE:
                Pair.append("""pair_coeff	 %i %i %s 0. %.2f  %.2f\n"""%(t1+1,num_particle["spb"]+n_i,precise,dist,cut))
        Pair.append("""pair_coeff	 %i %i %s 0. %.2f  %.2f\n"""%(num_particle["spb"],num_particle["spb"],precise,dist,cut))
        if DUPLICATE:
            Pair.append("""pair_coeff	 %i %i %s 0. %.2f  %.2f\n"""%(num_particle["spb"]+n_i,num_particle["spb"]+n_i,precise,dist,cut))
            Pair.append("""pair_coeff	 %i %i %s 0. %.2f  %.2f\n"""%(num_particle["spb"],num_particle["spb"]+n_i,precise,dist,cut))
        Pair.append("""pair_coeff	 %i %i gauss/cut  %.1f %.2f  %.2f \n"""%(num_particle["centro"],num_particle["spb"],-18000*2*radius/32., 0.3/realsigma,2*radius))
        if DUPLICATE:
            Pair.append("""pair_coeff	 %i %i gauss/cut  %.1f %.2f  %.2f \n"""%(num_particle["centro"]+n_i,num_particle["spb"]+n_i,-18000*2*radius/32., 0.3/realsigma,2*radius))
    
    g = open(REP+name,"w")
    g.write("".join(Bond)+"\n")
    g.write("".join(Pair)+"\n")
    if langle:
        g.write("".join(Angle))
    return spbond

def write_variable(ribosig,cellradius,maxi,etelo=5.0,ctelo=2.,name="variables",damp=2.0,diameter={},REP=""):
    f = open(REP+name,"w")    
    
    f.writelines("""
#########################################################
#initialisation variable
variable rad equal %.2f
variable mrad equal %.2f

#Physical parameter

variable sigDna equal %.2f
variable isigDna equal %.2f
variable sigrDNA equal %.2f
variable isigrDNA equal %.2f

variable frad equal %.2f    #Radius of the nucleus


#lj sigma parameter of the repulsif interation between the nucleus and DNA
variable eNorm equal 1.0
variable sigNorm equal ${sigDna}/2.
#Cutoff set to 1.12 this value; we keep only the repulsive part
variable sigNormCut equal ${sigNorm}*1.12


#lj sigma parameter of the repulsif interation between the nucleus and DNA
variable eRibo equal 1.0
variable sigRibo equal ${sigrDNA}/2.
#Cutoff set to 1.12 this value; we keep only the repulsive part
variable sigRiboCut equal ${sigRibo}*1.12

#Parameter of the attractive interaction between telomers and nucleus
variable etelo equal %.2f
variable sigtelo equal %.2f*${sigDna}
variable cuttelo equal %.2f*${sigDna}

#Damping
variable damp equal %.1f
"""%(maxi,maxi,diameter["1"],1./diameter["1"],diameter["3"],1./diameter["3"],cellradius,etelo,ctelo,3*ctelo,damp))
    f.close()

if __name__ == "__main__":
    import sys
########################################################
#Coarsed Model parameters
    #Fixed parameters
    compaction = 83.

    
    #Fixed length scale   
    sig=1 
    
    #Simu parameter:
    coarse = int(sig*float(sys.argv[2]))
    #Important parameters
    coarseribo=5000. # to generate 150 unit
    fractionrDNA = 1/6.      #Real volume is about twice    
    khun=1.
    cutoff=1.4
    
    khun = float(sys.argv[3])
    fractionrDNA = float(sys.argv[4])
    cutoff = float(sys.argv[5])
    etelo = float(sys.argv[6])
    ctelo = float(sys.argv[7])
    damp =  float(sys.argv[8])
    typecut = int (sys.argv[9])
    compaction = int (sys.argv[10])
    microtubule = float (sys.argv[11])
    cutlp=10
    REP =sys.argv[12] + "/"
    reducedv_factor = float(sys.argv[13])
    #secondary parameters
    

    #Deduced parameter    
    realsigma = coarse/83. / 1000. #in micrometer
    radius = sig/realsigma
    ribosig = 2*radius * (fractionrDNA/150)**(1/3.) 
    #radius = .95 * radius


    print ribosig , radius,150* (ribosig/2/radius)**3 * 6 
    print coarse,compaction
    print "density", 2400*coarse/5000 * compaction/83. /(5/6.*(2*radius)**3)
    #exit()

    cell = sys.argv[1]
    if cell == "poly":
        diameter={"1":sig,"3":sig,"2":sig}
        num_particle = {"telo":2,"ribo":3,"centro":4,"spb":5}
        SPB=False
        n=10
        angle_bond=False
        chrlen = [ 100*5000 for i in range(n)]
        chrcen = [50*5000 for i in range(n)]
        Ribopos=[]
        d = 2*0.056
        R = (n * 100/d) **(1/3.)/2.
        Extra = []

    if cell == "plasmo":

        diameter={"1":sig,"2":sig,"3":ribosig,"4":sig}
        num_particle = {"telo":2,"ribo":3,"centro":4,"spb":5}
        SPB=False
        angle_bond=True
        chrlen=[643292,
        947102,
        1060087,
        1204112,
        1343552,
        1377956,
        1350452,
        1323195,
        1541723,
        1694445,
        2035250,
        2271477,
        2747327,
        3291006]
        chrcen=[458021,
                446771,
                593691,
                648001,
                454541,
                477751,
                863741,
                299011,
                1241071,
                935681,
                830781,
                1281511,
                1167271,
                1070851]
    
        #generate pdb
        
        #Ribopos=[(1,473739,10000),(5,1289601,10000),(7,1083761,10000),(11,1923000,10000),(13,2700004,10000),(14,779377,10000)]
        #Ribopos=[(1,473739,125000),(5,1289601,125000),(7,1083761,125000),(11,1923000,125000),(13,2700004,125000),(14,779377,125000)]
        Extra=[]
        Ribopos=[(1,473739,10000),(5,1289601,375000),(7,1083761,375000),(11,1923000,10000),(13,2700004,10000),(14,779377,10000)]
    if cell == "yeastlike50":
        diameter={"1":sig,"2":sig,"3":ribosig,"4":sig}
        num_particle = {"telo":2,"ribo":3,"centro":4,"spb":5}
    
        SPB=False
        angle_bond=True
        ring=True
        Extra=[]
        Nch = 2500/50
        chrcen=[None for i in range(Nch)]
        
        chrlen=[50*5000 for i in range(Nch)]
        if ring:
            chrlen.append(16000)
            chrcen.append(None)
            chrlen.append(6000)
            chrcen.append(None)
            
        Ribopos=[(25,25*5000,750000)]
        
    if cell == "yeastlike100":
        diameter={"1":sig,"2":sig,"3":ribosig,"4":sig}
        num_particle = {"telo":2,"ribo":3,"centro":4,"spb":5}
    
        SPB=False
        angle_bond=True
        ring=True
        Extra=[]
        Nch = 2500/100
        chrcen=[None for i in range(Nch)]
        
        chrlen=[100*5000 for i in range(Nch)]
        if ring:
            chrlen.append(16000)
            chrcen.append(None)
            chrlen.append(6000)
            chrcen.append(None)
            
        Ribopos=[(12,50*5000,750000)]
        Ribopos=[]
    if cell == "yeast4only":
        diameter={"1":sig,"2":sig,"3":ribosig,"4":sig}
        num_particle = {"telo":2,"ribo":3,"centro":4,"spb":5}
    
        SPB=True
        if int(khun) == 0:
            angle_bond=False
        else:
            angle_bond=True
            
        ring=False
        Extra=[]
        chrlen=[1531933,
                20000]
        
        chrcen=[None,None]
       
        Ribopos=[(2,5000,5000)]
    if cell == "yeast":
        Mass = {"1":1,"2":1,"3":34,"4":1,"5":1}

        diameter={"1":sig,"2":sig,"3":ribosig,"4":sig}
        num_particle = {"telo":2,"ribo":3,"centro":4,"spb":5}
    
        SPB=True
        if khun == 0:
            angle_bond=False
        else:
            angle_bond=True
            
        ring=True
        Extra=[]
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
        if ring:
            chrlen.append(16000)
            chrcen.append(None)
            chrlen.append(6000)
            chrcen.append(None)
            
        Ribopos=[(12,450000,750000)]
        cutlp = 10
        if typecut == 1:
            rigidcut=6
            diameter["%i"%rigidcut] = sig
            Extra = [(4,161,10 ,rigidcut)]  #[4,"msd","ars425",161]
            chrcen[3] += 10 * 5000
        if typecut == 2:
            rigidcut=6
            diameter["%i"%rigidcut] = ribosig
            Extra = [(4,161,1 ,rigidcut)]  #[4,"msd","ars425",161]
            #chrcen[3] += 1 * 5000
        if typecut == 3:
            rigidcut=6
            diameter["%i"%rigidcut] = sig
            Extra = [(4,161,2 ,rigidcut)]  #[4,"msd","ars425",161]
        """
        if typecut == 4:
            rigidcut=6
            diameter["%i"%rigidcut] = 0.5*sig
            Extra = [(4,161,20 ,rigidcut)]
            chrcen[3] += 20 * 5000
            cutlp=1
        """    
        if typecut == 4:
            rigidcut=6
            diameter["%i"%rigidcut] = 0.1*sig
            Extra = [(4,161,20 ,rigidcut)]
            chrcen[3] += 100 * 5000
            cutlp=100
        if typecut == 5:
            cutlp=10
            rigidcut=6
            diameter["%i"%rigidcut] = sig
            Extra = [(4,161,20 ,rigidcut)]  #[4,"msd","ars425",161]
            chrcen[3] += 20 * 5000
        if typecut == 6:
            cutlp=100
            rigidcut=6
            diameter["%i"%rigidcut] = sig
            Extra = [(4,161,10 ,rigidcut)]  #[4,"msd","ars425",161]
            chrcen[3] += 10 * 5000
        if typecut == 7:
            cutlp=1
            rigidcut=6
            diameter["%i"%rigidcut] = sig
            Extra = [(4,161,10 ,rigidcut)]  #[4,"msd","ars425",161]
            chrcen[3] += 10 * 5000
            
        if typecut == 8:
            cutlp=10
            rigidcut=6
            diameter["%i"%rigidcut] = sig
            Extra = [(4,161,1 ,rigidcut)]  #[4,"msd","ars425",161]
            chrcen[3] += 1 * 5000
        if typecut == 9:
            cutlp=10
            rigidcut=6
            diameter["%i"%rigidcut] = sig
            Extra = [(4,161,2 ,rigidcut)]  #[4,"msd","ars425",161]
            chrcen[3] += 2 * 5000
        if typecut == 10:
            cutlp=1
            rigidcut=6
            diameter["%i"%rigidcut] = sig
            Extra = [(4,229,1 ,rigidcut)]  #[4,"msd","ars425",161]
            chrcen[3] += 1 * 5000
    if cell == "testsp":
        diameter={"1":sig,"2":sig,"3":ribosig,"4":sig}
        num_particle = {"telo":2,"ribo":3,"centro":4,"spb":5}
    
        SPB=False
        angle_bond=True
        ring=False
        Extra=[]
        chrlen=[50*5000,100*5000,0]
        #chrlen = [50*coarse]
        chrcen=[None,None,None]
        diameter["3"] = 1.1*numpy.sqrt(2)
        #num_particle["twic"] = 6
        Extra = [(1,50,25 ,3),(3,0,50,3)] 
        Ribopos=[]
            #chrcen[3] += 1 * 5000
    if cell == "testb":
        ribosig = 1.1*numpy.sqrt(2)

        diameter={"1":sig,"2":sig,"3":ribosig,"4":sig}
        Mass = {"1":1,"2":1,"3":2,"4":1,"5":1}
        num_particle = {"telo":2,"ribo":3,"centro":4,"spb":5}
    
        SPB=True
        angle_bond=True
        ring=True
        Extra=[]
        NN=10
        radius=10
        chrlen=[0*5000 for i in range(NN)]
        
        chrcen=[None for i in range(NN)]
        Ribopos=[]
        Extra=[(i+1,0,50,3)  for i in range(NN) ]
        
    if cell == "test":
        diameter={"1":sig,"2":sig,"3":ribosig,"4":sig}
        num_particle = {"telo":2,"ribo":3,"centro":4,"spb":5}
    
        SPB=True
        angle_bond=True
        ring=True
        Extra=[]
        NN=10
        radius=10
        chrlen=[100*5000 for i in range(NN)]
        
        chrcen=[None for i in range(NN)]
        Ribopos=[]
    if cell == "diplo":
        diameter={"1":sig,"2":sig,"3":ribosig,"4":sig}
        num_particle = {"telo":2,"ribo":3,"centro":4,"spb":5}
        #print ribosig

        SPB=True
        angle_bond=True
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
        948066,
        230218,
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
        555957,
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
        rigidcut=6
        Ribopos=[(12,450000,750000),(28,450000,750000)]
        diameter["%i"%rigidcut] = ribosig
        """
        Extra = [(5,23,10 ,rigidcut)]
        chrcen[4] += 10 * 5000
        """
        Extra = []
        radius = 1.5/realsigma
    
    chromosomes = [ int(ch/(coarse*compaction/83.)) if int(ch/(coarse*compaction/83.)) != 0 else 1  for ch in chrlen]
    print chromosomes
    centrolist  = [int(ch/(coarse*compaction/83.)) if ch != None else None for ch in chrcen]
    print centrolist
    import cPickle
    if Ribopos != []:
        Ribopos=[(a,int(b/(coarse*compaction/83.)),int(c/coarseribo)) for a,b,c in Ribopos]
    print Ribopos
    print sum(chromosomes)
    summary = [chromosomes,centrolist,Ribopos,Extra]
    f = open(REP+"summary",'w')
    cPickle.dump(summary,f)
    f.close()
    f = open(REP+"typecell","w")
    f.writelines("variable typecell index %s\n"%cell)
    f.close()
    if typecut == 3:
        f = open(REP+"scenari","w")
        f.writelines("""
dump 		1 all dcd 1000 dump_random.${typecell}.comp.dcd
run 10000000
#run 10000
group cut type %i
delete_bonds cut multi special 

undump 1

dump 		4 all dcd 1000 dump3.${typecell}.comp.dcd
dump 		2 all dcd 1 dumpshort.${typecell}.comp.dcd
dump 		3 all dcd 10 dumpint.${typecell}.comp.dcd

run 100000
undump 2
run 1000000
undump 3
run 100000000
unfix 1
unfix 2
write_restart 	restart.${simname}.dreiding1"""%rigidcut)
        f.close()
    elif typecut == 10:
        f = open(REP+"scenari","w")
        f.writelines("""
dump 		1 all dcd 1000 dump_random.${typecell}.comp.dcd
run 10000000
#run 10000
group cutted  type %i
fix wallcut cutted wall/region mySphere  lj93  ${etelo} ${sigtelo} ${cuttelo}



undump 1

dump 		4 all dcd 1000 dump3.${typecell}.comp.dcd
dump 		2 all dcd 1 dumpshort.${typecell}.comp.dcd
dump 		3 all dcd 10 dumpint.${typecell}.comp.dcd

run 100000
undump 2
run 1000000
undump 3
run 100000000
unfix 1
unfix 2
write_restart 	restart.${simname}.dreiding1"""%rigidcut)
        f.close()
    else:
        f = open(REP+"scenari","w")
        f.writelines("""
        
dump 		1 all dcd 10000 dumpevery10.${typecell}.comp.dcd
run 5000000

#run 100000000

undump 1
#dump 		2 all dcd 1 dumpshort.${typecell}.comp.dcd
#run 100000
#undump 2
#dump 		3 all dcd 10 dumpint.${typecell}.comp.dcd
#run 1000000
#write_restart 	restart.${simname}.dreiding1 
""")
        f.close()
    coord,maxi,spbond,softspbond = get_coord(chromosomes,centrolist,Ribopos,diameter,Mass,
                                  num_particle,cell=cell,angle_bond=angle_bond,SPB=SPB,khun=khun,Extra=Extra,cutlp=cutlp,microtubule=microtubule,REP=REP,reducedv_factor=reducedv_factor)
                                #  dists=0,boxs=73.04)
    #coord,maxi = duplicate(chromosomes,centrolist,Ribopos,diameter,num_particle,cell=cell,angle_bond=angle_bond,SPB=SPB,khun=khun,Extra=Extra,cutlp=cutlp)
    write_variable(ribosig,radius,1.4*maxi,etelo,ctelo,damp=damp,diameter=diameter,REP=REP)
    
    if SPB:
        with open(REP+'spbinteraction','w') as spb:
            spb.writelines(spbond)
        with open(REP+'softspbinteraction','w') as spb:
            spb.writelines(softspbond)
    else:
        with open(REP+'spbinteraction','w') as spb:
            spb.writelines("""#No spb""")
				