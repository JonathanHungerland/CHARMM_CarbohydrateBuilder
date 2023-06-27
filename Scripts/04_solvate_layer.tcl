if {![info exists water_padding]} {
    puts stderr "Defaulting \"water_padding\" to 15.0 Angström"
    set water_padding 15.0
}

proc water_solvate { layer_file water_padding } {
    package require solvate
    solvate $layer_file.psf $layer_file.pdb -t $water_padding -s WT -o solvated_layer
    return "solvated_layer"
}

proc simple_ionize_NACL { file salt_concentration } {
    global cation anion temperature
    mol new $file.psf
    mol addfile $file.pdb
    
    set charge_to_compensate [vecsum [[atomselect top "all"] get charge]]
    if { $charge_to_compensate > 0 } {
        set neutralize_cation 0
        set neutralize_anion  [expr round($charge_to_compensate)]
    } else {
        set neutralize_cation [expr -1*round($charge_to_compensate)]
        set neutralize_anion  0
    }
    puts "Ensuring neutralization requires putting $neutralize_cation cations and $neutralize_anion anions."
    
    set numwaters [[atomselect top "water and element O"] num]
    set saltdouble [calculate_nions $temperature $numwaters $salt_concentration]
    set saltinteger [expr {round($saltdouble)}]
    set ncation [expr $saltinteger + $neutralize_cation]
    set nanion  [expr $saltinteger + $neutralize_anion]
    set integer_concentration [calculate_concentration $temperature [expr $numwaters - $ncation - $nanion] $saltinteger]
    puts stderr "Requiring $saltdouble ions to meet the necessary salt concentration which was rounded to $saltinteger."
    puts stderr "Paying attention to deletion of water molecules by autoionize, effective concentration of $integer_concentration will result."

    package require autoionize
    if { $ncation == 0 && $nanion == 0} { 
    } elseif { $cation == "SOD" && $anion == "CLA" } {
       autoionize -psf $file.psf -pdb $file.pdb -nions [list [list $cation $ncation] [list $anion $nanion]] -o simple_ionized
    } else {
       puts "ERROR: Unsupported kinds of ions! Rework charge compensation protocol to calculate nions correctly!"
       exit
    }
    return "simple_ionized"
}

proc water_placer { layer_file nwaters } {
    #Solvate the polymer-ions layer with a certain number
    #of water molecules.
    if {$nwaters==0} {
        puts stderr "No solvation because nwaters=0."
        return $file
    }
    package require solvate
    solvate $layer_file.psf $layer_file.pdb -t 7 -s WT -o solvated_layer
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
