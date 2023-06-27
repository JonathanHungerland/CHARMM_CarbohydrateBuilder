if {![info exists cation]} {
    puts stderr "Defaulting \"cation" to \"SOD\"."
    set cation "SOD"
}
if {![info exists anion]} {
    puts stderr "Defaulting \"anion" to \"CLA\"."
    set cation "CLA"
}


proc ion_placer { file ncation nanion } {
    global cation anion
    set with_cation [load_ions $file $ncation "SOD"]
    set with_anions [load_ions $file $nanion "CLA"]
    mol new $with_anions.psf waitfor all
    mol addfile $with_anions.pdb waitfor all
    position_ions $cation
    position_ions $anion
    set all [atomselect top "all"]
    $all writepdb placed_ions.pdb
    $all writepsf placed_ions.psf
    return "placed_ions"
}

proc load_ions { file number type } {
    global fileslocation
    if {$number == 0} {return $file}
    puts stderr "Adding $number ions of type $type."
    package require topotools
    set midlist {}
    set mol [mol new $file.psf waitfor all]
    mol addfile $file.pdb
    lappend midlist $mol
    for {set i 0} {$i<$number} {incr i} {
        set mol [mol new $fileslocation/$type.psf waitfor all]
        mol addfile $fileslocation/$type.pdb
        lappend midlist $mol
    }
    set mol [::TopoTools::mergemols $midlist]
    animate write psf added_ion_${type}.psf $mol
    animate write pdb added_ion_${type}.pdb $mol
    mol delete all
    return "added_ion_${type}"
}

proc position_ions { type } {
    puts stderr "Picking charged polymer atoms to counter ions of type $type."
    set ionsel [atomselect top "resname $type"]
    if {[$ionsel num] == 0} {return}
    set ioncharge [lindex [$ionsel get charge] 0]
    set polys [atomselect top "not water and not ions"]
    set polindex_list [$polys get index]
    set polcharge_list [$polys get charge]
    puts -nonewline stderr  "Putting $type close to: "
    set start_i 0
    foreach iind [$ionsel get index] {
        for {set i $start_i} {[llength $polindex_list] != 0} {incr i} {
            set listi [expr $i % [llength $polindex_list]]
            set index [lindex $polindex_list $listi]
            set charge [lindex $polcharge_list $listi]
            if { [expr $charge * $ioncharge] > 0 } { continue }
            if { [expr rand()] < [expr abs($charge)] } {
                puts -nonewline stderr  "$index "
                set psel [atomselect top "index $index"]
                set ion [atomselect top "index $iind"]
                $ion moveto [measure center $psel]
                set good 0
                for {set tries 0} {$tries < 100} {incr tries} {
                    $ion moveby [list [expr rand()] [expr rand()] [expr rand()]]
                    if {[ion_check $iind $ioncharge]} {set good 1; break}
                }
                if {! $good} {continue}
                puts -nonewline stderr "ACCEPT "
                set polindex_list [lreplace $polindex_list $listi $listi]
                set polcharge_list [lreplace $polcharge_list $listi $listi]
                incr start_i 500
                $psel delete
                $ion delete
                break
            }
        }
    }
    puts stderr ".. Done."
}

proc ion_check { ionindex ioncharge } {
    set near_sel [atomselect top "(within 3.5 of index $ionindex) and (not index $ionindex)"]
    set close_sel [atomselect top "(within 2.0 of index $ionindex) and (not index $ionindex)"]
    if {[expr [vecsum [$near_sel get charge]]*$ioncharge] < 0 && [$close_sel num] == 0} {
        $near_sel delete
        $close_sel delete; 
        return 1
    }
    $near_sel delete
    $close_sel delete
    return 0
}
