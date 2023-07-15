package require topotools
mol new Example2.psf
mol addfile Example2.pdb

topo writelammpsdata Example2_initial.lammpsdata

#now to include the restart files
mol new Example2.psf
mol addfile equilibrated.restart.coor type namdbin
set all [atomselect top "all"]
package require pbctools
pbc readxst equilibrated.restart.xsc
pbc get
$all writepdb Example2_equilibrated.pdb

topo writelammpsdata Example2_equilibrated.lammpsdata
