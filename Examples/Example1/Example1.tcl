#Goal of this script: Construct two randomly oriented Methanochondroitin polymers
#inside a box of ionized water.

#Here are the settings:
set build_method "angry_baby" ;#will randomly throw the polymers around until they don't clash
set space_per_polymer "12.0" ;#the box in which angry_baby throws. it will be a cube of side length
                              #space_per_polymer*number_of_polymers. the latter is read automatically 
                              #once all polymers have been built
set temperature 300 ;#simulation temperature, only relevant to get salt concentration perfectly right.
set water_padding 20 ;#how much distance in angström there is to the edge of the water box.
set cation "SOD" ;#the cation for ionization
set anion "CLA" ;#the anion for ionization 
set salt_concentration "0.15" ;#salt concentration in mol/L. 0.00 is a valid entry if you don't want salts.

#import all scripts that you need. will also set default values for all parameters that have not been given yet.
source ../../Scripts/01_polymer_construction.tcl
source ../../Scripts/02_layer_construction.tcl
source ../../Scripts/03_ionize_layer.tcl
source ../../Scripts/04_solvate_layer.tcl
source ../../Scripts/05_finalize.tcl
set fileslocation [file normalize "../../TopologyData"]
cd tmp

set residues [list "GALP" "GLCA"  "GALN" ]
set links    [list     "14"    "13"    "13" ]
set repeats 2
#these lists will make build_polymer construct a polymer of the structure
#Protonated Galactosamine --(link from C1 carbon to C4 carbon)--> Glucoronic Acid --(C1 carbon to C3 carbon)--> Galactosamine --(C1 carbon to C3 carbon)--> #Protonated Galactosamine --(link from C1 carbon to C4 carbon)--> Glucoronic Acid --(C1 carbon to C3 carbon)--> Galactosamine
#valid residues are: GALN = Neutral Galactomasine                          \_in polymer environment probably 50% positive (pKa comparable to Lysine sidechain)
#                    GALP = Protonated Galactosamine (positively charged)  /
#                    GLCA = Glucoronic Acid (negatively charged)           \_probably mostly negative (pKa comparable to Glutamic Acid sidechain)
#                    GLCP = Protonated Glucoronic Acid (neutral)           /


set poly_list {}
for {set i 1} {$i <= 2} {incr i} {
   lappend poly_list [build_polymer $i $residues $links $repeats]
}
set layer [build_layer $poly_list]
set solvated [water_solvate $layer $water_padding]
set ionized [simple_ionize_NACL $solvated $salt_concentration]
set final [finalize $ionized]
file copy -force $final.psf ../Example1.psf
file copy -force $final.pdb ../Example1.pdb
file copy -force $final.pbc ../Example1.pbc

quit
