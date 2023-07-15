#Goal of this script: Construct a stacked structure of a polymer with some
#water molecules and some ions in between

#Here are the settings:
set build_method "bored_bricklayer" ;#will construct layers of parallel polymers
set brick_increment 1.0
#### POLY -> POLY -> POLY ####    ==> Columns
#### POLY -> POLY -> POLY ####  ||
#### POLY -> POLY -> POLY ####  \/ Rows
#### POLY -> POLY -> POLY ####    
set ncolumns 3 ;#number of columns
set nrows 8 ;#number of rows
set nstacks 8 ;#number of vertical stacks

set temperature 300 ;#simulation temperature, only relevant to get salt concentration perfectly right.
set cation "SOD" ;#the cation for ionization shall be sodium
set ncation 0
set anion "CLA" ;#the anion for ionization shall be chloride
set nanion 0
set nwaters [expr 50*$nrows*$ncolumns*$nstacks]

#import all scripts that you need. will also set default values for all parameters that have not been given yet.
source ../../Scripts/01_polymer_construction.tcl
source ../../Scripts/02_layer_construction.tcl
source ../../Scripts/03_ion_placer.tcl
source ../../Scripts/04_water_placer.tcl
source ../../Scripts/05_finalize.tcl
set fileslocation [file normalize "../../TopologyData"]
cd tmp

set residues [list "GALP" "GLCA"  "GALN" ]
set links    [list     "14"    "13"    "13" ]
set repeats 3


set poly_list {}
for {set i 1} {$i <= [expr $ncolumns*$nrows*$nstacks]} {incr i} {
   lappend poly_list [build_polymer $i $residues $links $repeats]
}
set layer [build_layer $poly_list]
set ionized [ion_placer $layer $ncation $nanion]
set solvated [water_placer $ionized $nwaters]
set final [finalize $solvated]
file copy -force $final.pbc ../Example2.pbc
set shifted_stacks [shift_stacks $final]
file copy -force $shifted_stacks.psf ../Example2.psf
file copy -force $shifted_stacks.pdb ../Example2.pdb


quit
