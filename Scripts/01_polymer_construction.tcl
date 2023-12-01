if {![info exists clash_tolerance_dist]} {
    puts stderr "Defaulting \"clash_tolerance_dist\" for polymer-internal clashes to 0.8 Angström."
    set clash_tolerance_dist 0.8
}
if {![info exists mindist]} {
    puts stderr "Defaulting \"mindist for next-residue-placement to 4.5 Angström."
    set mindist 4.5
}
if {![info exists maxdist]} {
    puts stderr "Defaulting \"maxdist\" for next-residue-placement to 6.5 Angström."
    set maxdist 6.5
}
if {[info exists persistence_length]} {
    puts stderr "persistence_length based residue placement enabled."
    puts stderr "Defaulting \"allowed_rotation\" for next-residue-placement to 40 degrees."
    puts stderr "And drawing from there with a gaussian distribution according to persistence_length $persistence_length"
} elseif {![info exists allowed_rotation]} {
    puts stderr "Defaulting \"allowed_rotation\" for next-residue-placement to 15 degrees."
    set allowed_rotation 15
}
if {![info exists max_residue_tries]} {
    puts stderr "Defaulting \"max_residue_tries\" for placing the next residue to 30."
    set max_residue_tries 50
}
if {![info exists bglc_protonation_prob]} {
    puts stderr "Defaulting probability \"bglc_protonation_prob\" for protonating residue BGLC to 0.0"
    set bglc_protonation_prob 0.0
}
if {![info exists polbuild_fromsel_text]} {
    puts stderr "Defaulting selection text \"polbuild_fromsel_text\" to \"name C3 O3 C4 O4\""
    set polbuild_fromsel_text "name C1 O1 HO1 H1"
}
if {![info exists polbuild_tosel_text]} {
    puts stderr "Defaulting selection text \"polbuild_tosel_text\" to \"name C1 O1 HO1 H1\""
    set polbuild_tosel_text "name C3 O3 C4 O4"
}
if {![info exists nrepeats]} {
    puts stderr "Defaulting number \"nrepreats\" of basic unit polymer repeats to 3."
}


proc residue_link_generator { residue_repeat link_repeat nrepeats } {
    #A possibly extendable function that spits out a list of
    #residues to construct the polymer and how to link them
    global bglc_protonation_prob all_residue_types
    
    set residuelist {}
    set linklist [list "skip"]
    set patchlist {}
    for {set i 0} {$i<$nrepeats} {incr i} {
        foreach residue $residue_repeat link $link_repeat {
            if { $residue == "GALN" || $residue == "GLCA" } {
                lappend residuelist $residue
                lappend patchlist "NONE"
            } elseif { $residue == "GALP" } {
                lappend residuelist "GALN"
                lappend patchlist "GALP"
            } elseif { $residue == "GLCP" } {
                lappend residuelist "GLCA"
                lappend patchlist "GLCP"
            } elseif { $residue == "N2P" } {
                lappend residuelist "GALN"
                lappend patchlist "N2P"
            } else {
                puts stderr "ERROR: Unrecognized residue type."
                exit
            }
            lappend linklist "${link}bb"
        }
    }
    set linklist [lreplace $linklist end end]

    lappend all_residue_types [lsort -u $residuelist]
    lappend all_residue_types [lsort -u $patchlist]
    return [list $residuelist $linklist $patchlist]
}

#MAIN FUNCTION
proc build_polymer { polnumber residuelist linklist nrepeats } {
    #Builds a polymer pdb and psf file based on given list of 
    #residues and how they shall be linked.

    global fileslocation max_residue_tries
    set generated [residue_link_generator $residuelist $linklist $nrepeats]
    set residueslist [lindex $generated 0]
    set linklist [lindex $generated 1]
    set patchlist [lindex $generated 2]
    
    puts stderr "Building polymer: $residueslist"
    puts stderr "With patches    : $patchlist"
    puts stderr "Linked by       : $linklist"

    #load in a bunch of residues as "building blocks"
    set lego_set [init_lego_set $fileslocation $residueslist $polnumber]
    mol new $lego_set
    place_first_residue 1 [lindex $residueslist 0]

    set tries_sum 0
    for {set r 1} {$r < [llength $residueslist]} {incr r} {
        set next_resid [expr $r + 1]
        set next_resname [lindex $residueslist $r]
        set tries [place_next_residue $next_resid $next_resname]
        if {$tries == 0} {
            puts stderr "No succes for placing residue after ${max_residue_tries} tries. Starting anew."
            return [build_polymer $polnumber]
        }
        set tries_sum [expr $tries_sum + $tries]
    }
    puts stderr "Polymer residues of segid P$polnumber were, on average, succesfully placed after [expr {double($tries_sum) / [llength $residueslist]}] tries."
    #write pdb function
    set all [atomselect top "all"]
    $all writepdb P$polnumber.pdb
    $all delete 

    set psfgen_file [polymer_psfgen P$polnumber $linklist $patchlist $polnumber]
    #psfgen function
    return $psfgen_file
}

proc polymer_psfgen { inputfile linklist patchlist polnumber } {
    global fileslocation
    package require psfgen
    resetpsf
    mol new $inputfile.pdb
    [atomselect top "all"] set segid P$polnumber

    topology $fileslocation/top_all36_carb_new.rtf

    topology alias GALN BGALNA
    topology alias GLCA BGLCA

    pdbalias atom GALN N2 N
    pdbalias atom GALN HN2 HN
    pdbalias atom GALN C7 C
    pdbalias atom GALN O7 O
    pdbalias atom GALN C8 CT
    pdbalias atom GALN H81 HT1
    pdbalias atom GALN H82 HT2
    pdbalias atom GALN H83 HT3

    pdbalias atom GLCA O6A O61
    pdbalias atom GLCA O6B O62

    segment P$polnumber {
        first NONE
        last NONE
        pdb $inputfile.pdb
    }
    coordpdb $inputfile.pdb P$polnumber

    #patch residues to get coorrect residue types (e.g. GALP instead of GALN)
    set r 0
    foreach patch $patchlist {
        incr r
        if {$patch=="NONE"} {continue}
        patch $patch P${polnumber}:$r
    }

    #link the residues
    set r 1
    foreach link $linklist {
        set r_before $r
        incr r
        if { $link == "skip" } {continue}
        puts "Doing link: ${link} P${polnumber}:$r_before P${polnumber}:$r"
        patch ${link} P${polnumber}:$r_before P${polnumber}:$r
    }
    puts "Links done!" 
    guesscoord
    regenerate angles dihedrals

    writepsf ${inputfile}_psfgen.psf
    writepdb ${inputfile}_psfgen.pdb
    delete_spurious_FEP_angles_dihedrals ${inputfile}_psfgen
    return "${inputfile}_psfgen"
}

proc delete_spurious_FEP_angles_dihedrals { inputfile } {
    #procedure to delete angles and dihedrals that regenerate angles dihedrals
    #has placed between vanishing/appearing atoms, which should not have interactions
    #with each other
    mol new $inputfile.psf
    mol addfile $inputfile.pdb
    package require topotools
    foreach angle [::TopoTools::topo getanglelist] {
       set angle_sel [atomselect top "index [lrange $angle 1 3]"]
       set atom_names [lsort -u [$angle_sel get name]]
       set A_found 0; set B_found 0
       foreach name $atom_names {
          if { [string match {A*} $name] } { set A_found 1 }
          if { [string match {B*} $name] } { set B_found 1 }
       }
       if { $A_found && $B_found } {
           puts "Deleting angle between indices [lrange $angle 1 3] with names $atom_names"
          ::TopoTools::topo delangle [lindex $angle 1] [lindex $angle 2] [lindex $angle 3]
       }
    }
    foreach dihedral [::TopoTools::topo getdihedrallist] {
       set dihedral_sel [atomselect top "index [lrange $dihedral 1 4]"]
       set atom_names [lsort -u [$dihedral_sel get name]]
       set A_found 0; set B_found 0
       foreach name $atom_names {
          if { [string match {A*} $name] } { set A_found 1 }
          if { [string match {B*} $name] } { set B_found 1 }    
       }
       if { $A_found && $B_found } {
          puts "Deleting dihedral between indices [lrange $dihedral 1 4] with names $atom_names"
          ::TopoTools::topo deldihedral [lindex $dihedral 1] [lindex $dihedral 2] [lindex $dihedral 3] [lindex $dihedral 4]
       }
    }  
    set all [atomselect top "all"]
    $all writepsf $inputfile.psf
    $all writepdb $inputfile.pdb
}

proc init_lego_set { pdblocation residueslist polnumber } {
    #Initialize a pdb file from which polymers are constructed from
    global polbuild_fromsel_text polbuild_tosel_text
    package require topotools
    set lego_pos {900 900 900}
    set sellist {}
    set i 0
    foreach resname $residueslist {
        mol new "${pdblocation}/${resname}.pdb" waitfor all
        set sel($i) [atomselect top "all"]
        $sel($i) moveby [vecinvert [measure center $sel($i)]]
        #align the molecule so that from_sel and to_sel align with the x axis
        set from_sel [atomselect top "$polbuild_fromsel_text"]
        set to_sel   [atomselect top "$polbuild_tosel_text"]
        set fv [measure center $from_sel]
        set tv [measure center $to_sel ]
        set dirvec [vecsub $tv $fv]
        set M [transvecinv $dirvec]
        $sel($i) move $M
        $from_sel delete; $to_sel delete

        #put the residue to its lego position
        $sel($i) moveby $lego_pos
        $sel($i) set resid 900
        $sel($i) set resname $resname
        lappend sellist $sel($i)
        set lego_pos [expr [lindex $lego_pos 0] + 10]
        set lego_pos [list $lego_pos $lego_pos $lego_pos]
        incr i
    }
    set mol [::TopoTools::selections2mol $sellist]
    animate write pdb lego_set_$polnumber.pdb $mol
    for {set j 0} {$j<$i} {incr j} {
        $sel($j) delete
    }
    mol delete all
    return "lego_set_$polnumber.pdb"
}

proc random_rotation_around_build_direction { selection } {
    set rand_degrees [expr rand()*360]
    $selection move [transaxis x $rand_degrees]
}

proc place_first_residue { resid resname } {
    #Choose the appropriate type of build-residue, lego-set-residues have resid 900
    #and have been moved by +900 +900 +900 or more.
    set legosel [atomselect top "resid 900 and resname $resname"]
    if { [$legosel num] == 0 } { puts stderr "ERROR: Can't find first residue to be built." }
    set ressel [atomselect top "same residue as (index [lindex [$legosel get index] 0])"]
    #Move the next residue to the center dedicated position
    $ressel moveby [vecinvert [measure center $ressel]]
    random_rotation_around_build_direction $ressel
    #Give the placed residue the correct resid
    $ressel set resid $resid
    $legosel delete; $ressel delete
}

proc place_next_residue { resid resname } {
    #Set the next polymer residue to a reasonable position
    global mindist maxdist allowed_rotation max_residue_tries
    global polbuild_fromsel_text polbuild_tosel_text
    global persistence_length

    #Determine a vector from the direction in which the former
    #residue points
    set last [expr $resid - 1]
    set lastsel [atomselect top "resid $last"]
    set lastcenter [measure center $lastsel]
    set from [atomselect top "resid $last and ($polbuild_fromsel_text)"]
    set to   [atomselect top "resid $last and ($polbuild_tosel_text)"]
    set tonext [vecnorm [vecsub [measure center $to] [measure center $from]]]

    #Choose the appropriate type of build-residue, lego-set-residues have resid 900
    #and have been moved by +900 +900 +900 or more.
    set legores [atomselect top "resid 900 and resname $resname"]
    if { [$legores num] == 0 } { puts stderr "ERROR: Can't find next residue to be built."}
    set ressel [atomselect top "same residue as index [lindex [$legores get index] 0]"]
    $ressel moveby [vecinvert [measure center $ressel]]
    random_rotation_around_build_direction $ressel
    $ressel set resid $resid

    set clash 1
    set try_count 0
    while {$clash} {
        incr try_count
        if {$try_count > ${max_residue_tries}} {return 0}
        #Next position is a vector at most 45 degrees rotated and between
        #mindist and maxdist far from the initial vector
        set scale [expr {$mindist + rand() * ($maxdist - $mindist)}]
        if {[info exists persistence_length]} {
            set next_build_direction [persistance_based_residue_placing $tonext $scale]
        } else {
            set next_build_direction [random_point_xdegrees_close $tonext $allowed_rotation]
        }
        set nextcenter [vecadd [vecscale $scale $next_build_direction] $lastcenter]

        #random rotation of the next residue about the axis of the next_build_direction
        $ressel move [transabout $next_build_direction [expr rand()*0.5] pi]
        $ressel moveby [vecinvert [measure center $ressel]]
        $ressel moveby $nextcenter
        set clash [internal_clash_test $resid]
    }
    $lastsel delete; $from delete; $to delete; $legores delete; $ressel delete;
    return $try_count
}

proc random_point_xdegrees_close { vector degree_dist }  {
    #Generate point on the unit sphere which is not more than
    #degree_dist rotational degrees away from the input vector
    #Highest allowed input is 45 degrees.

    set PI [expr acos(-1)]
    while {1} {
        #Generate point on unit sphere
        set u1 [expr rand()]; set u2 [expr rand()]
        set lambda [expr acos(2*$u1-1)-$PI/2]
        set phi [expr 2*$PI*$u2]
        set x [expr cos($lambda)*cos($phi)]
        set y [expr cos($lambda)*sin($phi)]
        set z [expr sin($lambda)]
        #Switch the sign to only consider close octants to the input vector. reduces necessary tries.
        #If e.g. x>1/sqrt(3) and the sign is different, then without a switch
        #the angle between the input-vector and the random vector would always
        #be larger than 45 degrees and only switching ensures this.
        set one_over_sqrt3 0.57733503
        if {[expr $x * [lindex $vector 0]]<0 && $x > $one_over_sqrt3} {set x [expr -$x]}
        if {[expr $y * [lindex $vector 1]]<0 && $x > $one_over_sqrt3} {set y [expr -$y]}
        if {[expr $z * [lindex $vector 2]]<0 && $x > $one_over_sqrt3} {set z [expr -$z]}
        set randvec [list $x $y $z]
        #Calculate angle between input and random vector
        set dot [vecdot $vector $randvec]
        #The random vector is already normalized here, otherwise the formula was different
        set angle [expr acos($dot)/[veclength $vector]]
        #Accept rotations smaller or equal to degree_dist
        if {[expr $angle*180/$PI] <= $degree_dist} {
            return $randvec
        }
    }
}

proc persistance_based_residue_placing { vector displacement } {
    #uses given information about the persistence length to generate a polymer
    #has a displacement appropriate to the expected persistence length
    #see: The Equilibrium Theory of Inhomogeneous Polymers by Glenn H. Fredrickson
    #page 47 eq. 2.74-2.75
    global persistence_length allowed_rotation
    set pstddev [expr $displacement / $persistence_length]
    while {1} {
        set try_vec [random_point_xdegrees_close $vector $allowed_rotation]
        set delta_u2 [vecdot [vecsub $vector $try_vec] [vecsub $vector $try_vec]]
        if { rand() < [get_gaussian $delta_u2 $pstddev] } {
            return $try_vec
        }
    }
}

proc get_gaussian { x stddev } {
    #generate point of a gaussian distribution with standard deviation stddev
    set PI [expr acos(-1)]
    set norm [expr 1.0/($stddev *sqrt(2*$PI))]
    set epart [expr exp(-1.0/2*($x/$stddev)^2)]
    return [expr $norm*$epart]
}

proc internal_clash_test { resid } {
    #Perform a clash test considering a given residue
    #Assume that a all residues that have already been constructed
    #have others_resid < next_resid
    global clash_tolerance_dist

    set others [atomselect top "(within $clash_tolerance_dist of resid $resid) and (not resid $resid)"]
    if { [$others num] } {
        #If the selection is not empty, there is a clash
        $others delete
        return 1
    }
    $others delete
    return 0
}



