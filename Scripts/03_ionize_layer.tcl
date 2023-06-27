package require autoionize
if {![info exists do_we_neutralize]} {
    puts stderr "Defaulting \"do_we_neutralize\" to \"yes\"."
    set do_we_neutralize "yes"
}


proc ionize_layer { file target_temp waters_per_monomer salt_type_list salt_concentration_list } {
    set nwaters [calculate_nwaters $file $waters_per_monomer]
    set salt_doubles_list {}; set salt_integers_list {}; set salt_integer_concentration {}
    foreach type $salt_type_list concentration $salt_concentration_list {
        set saltdouble [calculate_nions $target_temp $nwaters $concentration]
        lappend salt_doubles_list $saltdouble
        set saltinteger [expr {round($saltdouble)}]
        lappend salt_integers_list $saltinteger
        lappend salt_integer_concentration [calculate_concentration $target_temp $nwaters $saltinteger]
        if { $saltdouble > 0.0 && $saltinteger == 0 } {
            puts stderr "WARNING: Non-zero salt concentration was rounded to zero salt ions."
        }
    }
    puts stderr "Ion numbers calculated."
    puts stderr "List of desired salts:          $salt_type_list"
    puts stderr "List of desired concentrations: $salt_concentration_list"
    puts stderr "List of number of salt ions:    $salt_integers_list"
    puts stderr "List of final concentration based on number of salts:"
    puts stderr "                                $salt_integer_concentration"
    puts stderr "Placing ions. Might take a bit."
    set ionize_input_list {}
    set salts_to_do 0
    foreach salt_type $salt_type_list salt_num $salt_integers_list {
        if { $salt_num > 0 } { set salts_to_do 1 }
        lappend ionize_input_list [list $salt_type $salt_num]
    }
    if { $salts_to_do == 0 } { 
        puts stderr "No ionization is to do."
        return $file 
    }

    set ionize_this $file
    if { ! [is_solvated $file] } {
        set ionize_this [pre_solvation $file]
        autoionize -psf $ionize_this.psf -pdb $ionize_this.pdb -nions $ionize_input_list -o pre_ionized_layer
        mol new pre_ionized_layer.psf
        mol addfile pre_ionized_layer.pdb
        set notwat [atomselect top "not water"]
        $notwat writepdb ionized_layer.pdb
        $notwat writepsf ionized_layer.psf
    } else {
        autoionize -psf $ionize_this.psf -pdb $ionize_this.pdb -nions $ionize_input_list -o ionized_layer
    }
    return "ionized_layer"
}

proc salt_neutralization { file cation anion } {
    global do_we_neutralize
    if {$do_we_neutralize == "no"} {return $file}
    set ionize_this $file
    if { ! [is_solvated $file] } {
        set ionize_this [pre_solvation $file]
        if { $cation == "H3O" && $anion == "OH" } {
            return [water_neutralization $ionize_this]
        }
        autoionize -psf $ionize_this.psf -pdb $ionize_this.pdb -neutralize -cation $cation -anion $anion -o pre_neutralized_layer
        mol new pre_neutralized_layer.psf
        mol addfile pre_neutralized_layer.pdb
        set notwat [atomselect top "not water"]
        $notwat writepdb neutralized_layer.pdb
        $notwat writepsf neutralized_layer.psf
    } else {
        if { $cation == "H3O" && $anion == "OH" } {
            return [water_neutralization $ionize_this]
        }
        autoionize -psf $ionize_this.psf -pdb $ionize_this.pdb -neutralize -cation $cation -anion $anion -o neutralized_layer
    }
    return "neutralized_layer"
}

proc water_neutralization { file } {
    if {$do_we_neutralize == "no"} {return $file}
    mol new $file.psf
    mol addfile $file.pdb
    set all [atomselect top "all"]
    set charge_int [expr round([vecsum [$all get charge]])]
    if { $charge_int == 0 } {
        return $file
    } elseif { $charge_int > 0 } {
        set neutralizer "OH"
    } elseif { $charge_int < 0 } {
        set neutralizer "H3O"
    } else {
        exit 4 "ERROR: Could not determine neutralizing water for charge $charge_int"
    }
    $all delete
    set water_oxys [atomselect top "water and element O"]
    set water_oxys_inds [$water_oxys get index]
    set numoxys [llength $water_oxys_inds]
    set seglist {}
    set reslist {}
    for {set c 0} {$c < $charge_int} {incr c} {
        set randoxy [lindex $water_oxys_inds [expr round(rand()*$numoxys)]]
        set sel [atomselect top "same residue as index $randoxy"]
        lappend seglist [$sel get segid]
        lappend reslist [$sel get resid]
        $water_oxys delete
        $sel delete
    }
    package require psfgen
    resetpsf
    readpsf $file.psf pdb $file.pdb
    topology $fileslocation/toppar_waters.str
    foreach seg $seglist res $reslist {
        patch P$neutralizer $seg:$res
    }
    writepsf ${neutralized}.psf
    writepdb ${neutralized}.pdb
    return "neutralized"
}
proc calculcate_water_density { temperature } {
    #Calculate the water density in mol/L.
    #From http://cecs.wright.edu/people/faculty/sthomas/htappendix01.pdf
    #density of water at temperatures given in degrees celcius
    set celsius_of_density { 0.01 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 }
    set density_water_kgperm3 { 999.8 999.9 999.7 999.1 998.0 997.0 996.0 994.0 992.1 990.1 988.1 985.2 983.3 980.4 977.5 974.7 971.8 968.1 965.3 961.5 957.9 }
    #conversion:
    set kelvin_of_density {}
    foreach T $celsius_of_density {lappend kelvin_of_density [expr $T + 273.15]}
    if {$temperature > [lindex $kelvin_of_density end]} {
        exit -666 "Desired setup temperature is above calculation range for water density. Can only go up to (100+273.15)K."
    }
    set molarmass_water 18.01528 
    #in g/mol 
    set density_water_molperL {}
    foreach d $density_water_kgperm3 {lappend density_water_molperL [expr $d / $molarmass_water]}

    #get correct density of water by linear interpolation of data points
    set water_density 0
    for {set Tindex 1} {$Tindex < [llength $kelvin_of_density]} {incr Tindex} {
        if { [lindex $kelvin_of_density $Tindex] > $temperature } {
            set leftTemp [lindex $kelvin_of_density [expr $Tindex - 1]]
            set rightTemp [lindex $kelvin_of_density $Tindex]
            set tempdist [expr $rightTemp - $leftTemp]
            set leftDens [lindex $density_water_molperL [expr $Tindex - 1]]
            set rightDens [lindex $density_water_molperL $Tindex]
            set densdist [expr $rightDens - $leftDens]
            set densSlope [expr $densdist / $tempdist]
            set water_density [expr $leftDens + $densSlope * ($temperature-$leftDens)]
            break
        }
    }
    return $water_density
}

proc calculate_nions { temperature nwaters salt_concentration } {
    #Calculate the number of ions based on salt concentrations in the environment
    #of a polymer matrix given in units of mol/L. The conversion is based on the density
    #of water at the given temperature. The same amount of ions per water molecules shall
    #be achieved.
    set water_density [calculcate_water_density $temperature]
    #we require [mol/L(SALT)]/[mol/L(WATER)]=Nsalt/Nwater --> Nions = Nwater * [mol/L(SALT)]/[mol/L(WATER)]
    return [expr {$nwaters * $salt_concentration / $water_density}]
}

proc calculate_concentration { temperature nwaters number_of_molecules } {
    #Calculate the concentration in mol/L of a molecule when a certain number of molecules
    #exists in the current simulation cell. The concentration is based assuming that the
    #concentration in the present water is the same as in the surrounding bulk water
    if { $nwaters == 0 } { return "!no_waters!" }
    set water_density [calculcate_water_density $temperature]
    return [expr {double($number_of_molecules) / $nwaters * $water_density}]
}

proc pre_solvation { file } {
    package require solvate
    solvate $file.psf $file.pdb -t 20 -s WT -o pre_solvation
    return "pre_solvation"
}

proc is_solvated { file } {
    set solvated 0
    mol new $file.psf
    mol addfile $file.pdb
    set wats [atomselect top "water"]
    if { [$wats num ] > 0 } { set solvated 1 }
    $wats delete
    return $solvated
}
