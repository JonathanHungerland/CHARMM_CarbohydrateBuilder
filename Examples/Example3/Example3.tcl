#Goal of this script: Construct a stacked structure of a polymer with some
#water molecules and some ions in between

#Here are the settings:
set build_method "bored_bricklayer" ;#will construct layers of parallel polymers
#### POLY -> POLY -> POLY ####    ==> Columns
#### POLY -> POLY -> POLY ####  ||
#### POLY -> POLY -> POLY ####  \/ Rows
#### POLY -> POLY -> POLY ####    
set ncolums 2 ;#number of columns
set nrows 5 ;#number of rows
set nstacks 5 ;#number of vertical stacks

set temperature 300 ;#simulation temperature, only relevant to get salt concentration perfectly right.
set cation "SOD" ;#the cation for ionization shall be sodium
set ncation 5
set anion "CLA" ;#the anion for ionization shall be chloride
set nanion 5
set nwaters 50

#import all scripts that you need. will also set default values for all parameters that have not been given yet.
source ../../Scripts/01_polymer_construction.tcl
source ../../Scripts/02_layer_construction.tcl
source ../../Scripts/03_ion_placer.tcl
source ../../Scripts/04_water_placer.tcl
source ../../Scripts/05_finalize.tcl
set fileslocation [file normalize "../../TopologyData"]
cd tmp

set residues [list "GLCA" "GALN"  "GALN" ]
set links    [list     "13"    "13"    "14" ]
set repeats 3


set poly_list {}
for {set i 1} {$i <= [expr $ncolumns*$nrows*$nstacks]} {incr i} {
   #each polymer will be initiated as "lego_set" of merged pdb containing all needed monomers.
   #afterwards, the monomers are placed on after another following the desired order.
   #build_polymer then returns the basename of the filelocation as a string.
   lappend poly_list [build_polymer $i $residues $links $repeats]
}
#build_layer 
set layer [build_layer $poly_list]
set ionized [ion_placer $layer $ncation $nanion]
set solvated [water_placer $ionized $nwaters]
set final [finalize $solvated]
file copy -force $final.psf ../Example3.psf
file copy -force $final.pdb ../Example3.pdb
file copy -force $final.pbc ../Example3.pbc

quit
