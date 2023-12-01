package require topotools
mol new Example2.psf
mol addfile Example2.pdb
topo writegmxtop Example2.itp [list "par_all36_carb.prm"]
