# CHARMM_CarbohydrateBuilder
This is a collection of scripts to construct randomly orientied and structured systems
of carbohydrate polymer structures. It produces pdb/psf files that are needed for Molecular 
Dynamics simulations. Currently, it is only implemented for constructing the Methanochondroitin
polymer but it can also be used for other carbohydrates after indlucing the necessary
source files and some alterations to 01_polymer_generation.tcl .

All scripts are supposed to be run using VMD that can be downloaded for free from  
https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD  
and then they should be executed as  
vmd -e Example1.tcl
