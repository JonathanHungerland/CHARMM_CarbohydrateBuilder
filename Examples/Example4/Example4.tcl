#import all scripts that you need. will also set default values for all parameters that have not been given yet.
source ../../Scripts/01_polymer_construction.tcl
source ../../Scripts/02_layer_construction.tcl
source ../../Scripts/03_ionize_layer.tcl
source ../../Scripts/04_solvate_layer.tcl
source ../../Scripts/05_finalize.tcl
set fileslocation [file normalize "../../TopologyData"]
cd tmp

set residues [list "N2P" ]
set links    [list        ]
set repeats 1

set poly [build_polymer 1 $residues $links $repeats]

package require solvate
solvate $poly.psf $poly.pdb -t 20 -o ../GALP_solvated
write_pbc ../GALP_solvated

package require autoionize
autoionize -sc 1.20 -psf ../GALP_solvated.psf -pdb ../GALP_solvated.pdb -o ../GALP_ionized
write_pbc ../GALP_ionized

quit
