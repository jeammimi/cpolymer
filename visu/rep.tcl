# VMD for LINUXAMD64, version 1.9.1 (February 1, 2012)
# Log file '/home/jarbona/lamp_pol/rep.tcl', created by user jarbona
#mol delrep 0 0
set pdbname [lindex $argv 0]
set dcdname [lindex $argv 1]
set confname [lindex $argv 2]
mol load pdb $pdbname dcd $dcdname 
topo clearbonds 0
topo readlammpsdata $confname molecular
set lc [topo getbondlist -molid 1]
topo setbondlist $lc -molid 0
mol delete 1
mol delrep 0 0
#mol new  $pdbname  autobonds off 
#mol addfile $dcdname
mol color Chain
mol representation CPK 1.500000 0.100000 10.000000 10.000000
mol selection all
mol material Opaque
mol addrep 0
mol color Chain
mol representation  CPK 4.900000 0.300000 10.000000 10.000000
mol selection name ribo
mol material Opaque
mol addrep 0
mol color Name
mol representation VDW .900000 12.000000
mol selection name telo
mol material Opaque
mol addrep 0
mol color Name
mol representation VDW .900000 12.000000
mol selection name cent or name spbb
mol material Opaque
mol addrep 0
mol color Name
mol representation VDW 0.900000 12.000000
mol selection name rcut
mol material Steel
mol addrep 0
display resetview
# VMD for LINUXAMD64, version 1.9.1 (February 1, 2012)
# end of log file.
