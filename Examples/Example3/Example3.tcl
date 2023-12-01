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
set ncation 0
set anion "CLA" ;#the anion for ionization shall be chloride
set nanion 0
set nwaters 0

#import all scripts that you need. will also set default values for all parameters that have not been given yet.
source ../../Scripts/01_polymer_construction.tcl
source ../../Scripts/02_layer_construction.tcl
source ../../Scripts/03_ion_placer.tcl
source ../../Scripts/04_water_placer.tcl
source ../../Scripts/05_finalize.tcl
set fileslocation [file normalize "../../TopologyData"]
cd tmp

set residues [list "GLCP" "GALN"  "GALN" ]
set links    [list     "13"    "13"    "14" ]
set repeats 3


set poly_list {}
for {set i 1} {$i <= [expr $ncolumns*$nrows*$nstacks]} {incr i} {
   #each polymer will be initiated as "lego_set" of merged pdb containing all needed monomers.
   #afterwards, the monomers are placed on after another following the desired order.
   #build_polymer then returns the basename of the filelocation as a string.
   lappend poly_list [build_polymer $i $residues $links $repeats]
}
#build_layer merges all polymers into one file and the call the necessary build_method.
#in the case of "bored_bricklayer", the +x, -y etc. in the stderr indicate the direction in which it was build.
#+x and -x alternate because otherwise we will get a staircase-like structure but we actually want something cube-like instead.
#starting from the current position, the polymer is pushed in the build-direction until there are no clashes and a litte "wiggle"
#induces a little more randomness into the process.
set layer [build_layer $poly_list]
set final [finalize $layer]
file copy -force $final.psf ../Example3.psf
file copy -force $final.pdb ../Example3.pdb
file copy -force $final.pbc ../Example3.pbc

quit
