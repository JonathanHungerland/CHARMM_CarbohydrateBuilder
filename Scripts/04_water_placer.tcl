proc water_placer { layer_file nwaters } {
    #Solvate the polymer-ions layer with a certain number
    #of water molecules.
    if {$nwaters==0} {
        puts stderr "No solvation because nwaters=0."
        return $file
    }
    package require solvate
    solvate $layer_file.psf $layer_file.pdb -t 1 -s WT -o solvated_layer
    mol new solvated_layer.psf
    mol addfile solvated_layer.pdb
    set watsel [atomselect top "water and element O"]
    set watsel_list [$watsel get index]
    $watsel delete

    #choose nwaters water molecules to keep.
    set chosenwaters {}
    for {set i 0} {$i < $nwaters} {incr i} {
        set index [lindex $watsel_list [expr {int(rand() * [llength $watsel_list])}]]
        if { [lsearch $chosenwaters $index] != -1 } {incr i -1; continue}
        lappend chosenwaters $index
    }
    set chosenwaters [flatten $chosenwaters]
    set chosen_water_and_layer [atomselect top "not water or (same residue as index $chosenwaters)"]
    $chosen_water_and_layer writepdb solvated_layer.pdb
    $chosen_water_and_layer writepsf solvated_layer.psf
    $chosen_water_and_layer delete
    return "solvated_layer"
}

proc calculate_nwaters { layer_file waters_per_monomer } {
    #Calculate the number of water molecules that should be present in the
    #system based on the expected hydration shell.
    mol new $layer_file.psf waitfor all
    mol addfile $layer_file.pdb waitfor all
    set polymers [atomselect top "not water and not ions"]
    set monomer_num 0
    foreach segid [lsort -u [$polymers get segid]] {
        set sel [atomselect top "segid $segid"]
        incr monomer_num [llength [lsort -u [$sel get resid]]]
        $sel delete
    }
    set nwaters [expr $waters_per_monomer * $monomer_num]
    return $nwaters
}

proc flatten {lst} {
    while 1 {
        set newlst [concat {*}$lst]
        if {$newlst == $lst} { break }
        set lst $newlst
    }
    return $newlst
}
