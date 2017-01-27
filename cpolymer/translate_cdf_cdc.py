# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 15:20:36 2013

@author: jarbona
"""
from .sortn import sort_nicely
import string
    
def write_single_pdb(name,chlist,namelist=[],origine=[],extraatom="",allzero=False):
    f = open(name,'w')
    alternate=" "
    code_insertion_residue=" "
    occupancy=1
    temp=0.0
    atom="ATOM  "
    Sum=0
    maxi = 0
    for n,(x,y,z,col) in enumerate(chlist):
        #print len(x)
        l = len(x)
        for i in range(l):
            
            number=i + Sum + 1
            name="bead"
            residue_name="bea"
            chain_id=string.ascii_letters[(n + 1) % len(string.ascii_letters)]
            residue_number=i
            X,Y,Z = x[i],y[i],z[i]
            if origine != []:
                X -= origine[n][0][i]
                Y -= origine[n][1][i]
                Z -= origine[n][2][i]
            if allzero:
                X,Y,Z=0,0,0
            maxi = max(X,Y,Z)
            seg_id="    "
            element_symbole=" "
            charge=" "
        #    print s+s+s+s+s+s+s+s

            traduction = {"1":"bead","2":"telo","3":"ribo","4":"cent","5":"spbb","6":"rcut"}
            name = traduction["%i" % namelist[n][i]]
            

            f.write("%6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"%(atom,number, name, alternate,residue_name,
        			 chain_id,residue_number,code_insertion_residue,X,Y,Z,occupancy,temp,
        				seg_id,element_symbole,charge))
        Sum += l
    if  extraatom != "":
        f.write(extraatom)
    Sum=1
    if maxi > 10000 or maxi < 0.0001:
        print(("Scale pb" , maxi)) 
        
    f.write("\n")
    for n,(x,y,z,col) in enumerate(chlist):
        l = len(x)
        for i in range(l - 1):
            number=i + Sum
            f.write("CONECT%5i%5i\n"%(number,number + 1))
        Sum += l
        
    f.close()
    
def write_single_xyz(name,chlist,origine=[]):
    fich = open(name,'w')
    f=["prout\n"]
    alternate=" "
    code_insertion_residue=" "
    occupancy=1
    temp=0.0
    atom="ATOM  "
    Sum=0
    maxi = 0
    for n,(x,y,z,col) in enumerate(chlist):
        l = len(x)
        #print len(x),len(y),len(z)
        for i in range(l):
            
            number=i + Sum + 1


            X,Y,Z = x[i],y[i],z[i]
            #print X,Y,Z
            if origine != []:
                X -= origine[n][0][i]
                Y -= origine[n][1][i]
                Z -= origine[n][2][i]


                
            f.append("C %30.10f %30.10f %30.10f\n"%(X,Y,Z))
            
        Sum += l
    f.insert(0,"%i\n"%Sum)

    fich.write("".join(f))
    fich.close()

def write_psf(name,chlist):
        f = open(name,'w')
        n = 0
        Sum = 0
        f.write("PSF CMAP\n")
        f.write("    1237 !NBOND: bonds\n")
        s=""
        for nc,c in enumerate(chlist):
           l = len(c[0])
           for i in range(l - 1):
               number=i + Sum + 1
               s +="%8i%8i"%(number,number + 1)
               n += 1
               if n == 4 or (nc == len(chlist) - 1 and i == l-2):
                   f.write(s+"\n")
                   n=0
                   s=""
           Sum += l
        f.close()
                   
               
                
if __name__ == "__main__":
        
    import glob,os,string
    from readm import read_name, return_ribo, return_centro, return_ch_info
    import sys
    S=1
    if len(sys.argv) == 3:
        S = int(sys.argv[2])
    
    root = "/home/jarbona/simunoyau/" + sys.argv[1] + "/"
    #root='/home/jarbona//cluster/faster-closing-applying-force-telomere-when-close-plasmo-bigger-ribo/'
    #root = 'whithout-microtubule-faster-closing-plasmo/'
    REDO=False
    SKIPalreadydone=False
    CLEAN=False
    center=False
    only_init=False
    all_zero=True
    
    
    rep = '%s/traj%i/'%(root,S)
    centro = 'parm6e-08_2.5e-08_11_2e-06_2e-07_9000000.%i'%S
    ribo = 'InitCoord6e-08_2.5e-08_11_2e-06_2e-07_9000000.%i'%S
    lc,tubule = return_ch_info(root+ribo)
    NC = len(lc)
    
    if os.path.exists(rep+'/scale'):
        f = open(rep+'/scale','r')
        scale = 1/float(f.readline())
        #scale = 1e7
        f.close()
    else:
        scale = 1e6
    if os.path.exists(rep+'/traj.dcd') and (SKIPalreadydone) and not(only_init):
        exit()

    
        
    test = glob.glob(rep+'/step*')#['step20000_30000']
    sort_nicely(test)
    to_concat = []
    for s,d in enumerate(test):
    
        d = d.split('/')[-1]
        where = rep+"/"+d+'/'
        ls = glob.glob(where+'traj*.gz')
        g_step = int(d.split('_')[0][4:])
        write_to = '%s/file_all%010i.dcd'%(where,g_step)
        print((where, s * 1.0 / len(test)))
        if os.path.exists(write_to) and (not REDO) and not(only_init):
            to_concat.append(write_to)
            continue
        origine=[]
       
        
    
        if len(ls) == 100 or (len(ls) == 99 and s == 0):
            sort_nicely(ls)

            for p,f in enumerate(ls):
                if s == 0 and p == 0 and center:
                    origine = read_name(f,rescale=scale)
                    print(scale)
                #try:
                ch = read_name(f,rescale=scale)
                """
                except:
                    print "Error reading" ,f
                    continue
                """
                step = int(f.split('.')[1][2:])
                pdb_file = rep+"/"+d+'/'+'%010i'%step+'.pdb'
                xyz_file = rep+"/"+d+'/'+'%010i'%step+'.xyz'
                #write_single_pdb(pdb_file,ch,centro_l,ribo_l,origine)
                #write_single_pdb(pdb_file,ch,centro_l,ribo_l,origine)
                write_single_xyz(xyz_file,ch,origine)
#                try:
#                    write_single_xyz(xyz_file,ch,origine)
#                except:
#                    print "Error reading" ,f
                if s == 0 and p == 0:
                    print("snssssnssn")
                    centro_l = return_centro(root+centro,nc=NC)
                    ribo_l = return_ribo(root+ribo,nc=NC)
                    print(ribo_l)
                    
                    write_single_pdb(pdb_file,ch,centro_l,ribo_l,origine,allzero=all_zero)
                    os.system('cp %s %s/init_conf.pdb'%(pdb_file,rep))
                
            #os.system('catdcd -o %s -pdb %s/*.pdb > /dev/null'%(write_to,where))
                if only_init:
                    break
            if only_init:
                    break
            os.system('catdcd -o %s -xyz %s/*.xyz > /dev/null'%(write_to,where))
            os.system('rm %s/*.xyz'%(where))
            to_concat.append(write_to)
    os.system('catdcd -o %s/traj.dcd  %s'%(rep," ".join(to_concat)))
    if CLEAN:
        os.system('rm %s'%(" ".join(to_concat)))
        
        