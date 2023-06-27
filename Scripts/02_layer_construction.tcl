if {![info exists space_per_polymer]} {
    puts stderr "Defaulting \"space_per_polymer\" to 10^3 Angström."
    set space_per_polymer 10
}
if {![info exists build_method]} {
    puts stderr "Defaulting \"build_method\" to angry_baby"
    set build_method angry_baby
}
puts stderr "Using build_method \"$build_method\""
if { $build_method == "angry_baby" } {
    if {![info exists maximum_layer_throws]} {
        puts stderr "Defaulting \"maximum_layer_throws\" for polymer layer to 100."
        set maximum_layer_throws 100
    }
    if {![info exists maximum_polymer_throws]} {
        puts stderr "Defaulting \"maximum_polymer_throws\" for single polymers to 100."
        set maximum_polymer_throws 1000
    }
    if {![info exists external_clash_tolerance]} {
        puts stderr "Defaulting \"external_clash_tolerance\" for clashes between polymers to 5.0 Angström."
        set external_clash_tolerance 5.0
    }
}
if { $build_method == "bored_bricklayer" } {
    if {![info exists ncolumns]} {
        puts stderr "Defaulting \"ncolumns\" to 2."
        set ncolumns 2
    }
    if {![info exists nrows]} {
        puts stderr "Defaulting \"nrows\" to 7."
        set nrows 7
    }
    if {![info exists nstacks]} {
        puts stderr "Defaulting \"nstacks\" to 3."
        set nstacks 3
    }
    if {![info exists brick_increment]} {
        puts stderr "Defaulting \"brick_increment\" to 0.5."
        set brick_increment 0.5
    }
}

proc minimize_polymers { layer_file namdexe } {
    global env fileslocation scriptslocation
    mol new $layer_file.psf waitfor all
    mol addfile $layer_file.pdb waitfor all
    set all [atomselect top "all"]
    set minmax [measure minmax $all]
    set pbcvec [vecadd [vecsub [lindex $minmax 1] [lindex $minmax 0]] {5.0 5.0 5.0}]
    set xlength [expr round([lindex $pbcvec 0])]; set ylength [expr round([lindex $pbcvec 1])];
    set zlength [expr round([lindex $pbcvec 2])];
    $all delete
    puts stderr "Minimizing polymer coordinates with NAMD2..."
    set namdexit [exec $env(SHELL) -c "$scriptslocation/03_minimize_layer.sh $namdexe $layer_file $fileslocation $xlength $ylength $zlength"]
    puts stderr "NAMD minimization of polymers exited with status $namdexit."
    mol delete top
    mol new $layer_file.psf waitfor all
    mol addfile "polymers_minimization/min_layer.dcd" waitfor all
    animate goto end
    set all [atomselect top "all"]
    $all writepdb ${layer_file}_minimized.pdb
    $all writepsf ${layer_file}_minimized.psf
    $all delete
    return "${layer_file}_minimized"
}

proc build_layer { polymers_list } {
    global build_method
    puts stderr "Building polymer layer."
    set merged_file [merge_polymers $polymers_list ]
    puts stderr "All polymers merged. Calling build method ${build_method}"
    set layer [$build_method $merged_file]
}

proc bored_bricklayer { file } {
    global npolymers
    global nrows ncolumns nstacks
    #stacks polymers next to and on top of one another
    mol new ${file}.psf waitfor all
    mol addfile ${file}.pdb waitfor all
    set npolymers [llength [lsort -u [[atomselect top "all"] get segid]]]
    #align all polymers in the x-y-plane pointing in x-direction
    for {set p 1} {$p<=$npolymers} {incr p} {
        set poly [atomselect top "segid P$p"]
        $poly moveby [vecinvert [measure center $poly]]
        set resids [lsort -dictionary -u [$poly get resid]]
        set start [lrange $resids 0 2]
        set end [lrange $resids end-2 end]
        set start_sel [atomselect top "segid P$p and resid $start"]
        set end_sel [atomselect top "segid P$p and resid $end"]
        set sv [measure center $start_sel]
        set ev [measure center $end_sel ]
        set dirvec [vecsub $sv $ev]
        set M [transvecinv $dirvec]
        $poly move $M
        #move them upwards to not clash later
        $poly moveby [list 0 0 [expr 300+$p*5]]
        $poly delete; $start_sel delete; $end_sel delete
    }

    set pID 1
    for {set stack 0} {$stack < $nstacks} {incr stack} {
        set pID [brick_stack $stack $pID]
    }
    set all [atomselect top "all"]
    $all writepdb bricklayer_build.pdb
    $all writepsf bricklayer_build.psf
    return "bricklayer_build"
}

proc brick_stack { stack pID } {
    puts stderr "Laying stack number $stack."
    global npolymers nrows ncolumns
    if {$pID > $npolymers} {return}
    if {$stack == 0} {
        set poly [atomselect top "segid P$pID"]
        $poly moveby [vecinvert [measure center $poly]]
        incr pID
    } else {
        set pID [lay_brick $pID "+z"]
    }
    for {set row 0} {$row < $nrows} {incr row} {
        set pID [brick_row $stack $row $pID]
    }
    xy_adjust_stack $stack
    return $pID
}

proc xy_adjust_stack { stack } {
    set stacksel [atomselect top "[get_stack $stack]"]
    if {$stack == 0} {
        $stacksel moveby [vecinvert [measure center $stacksel]]
    } else {
        set stack_center [measure center $stacksel]
        set prevsel [atomselect top "[get_stack [expr $stack - 1]]"]
        set prev_center [measure center $prevsel]
        set relative [vecsub $stack_center $prev_center]
        set move [vecinvert [vecmul {1 1 0} $relative]]
        $stacksel moveby $move
        $prevsel delete
    }
    $stacksel delete
    return
}

proc get_stack { stack } {
    #gives the atomselection for a given stack
    global nrows ncolumns
    set stack_string "segid "
    for {set s [expr 1 + $stack*$nrows*$ncolumns]} {$s<[expr 1+($stack+1)*$nrows*$ncolumns]} {incr s} {
        append stack_string "P$s "
    }
    return $stack_string
}

proc brick_row { stack row pID } {
    puts stderr "Laying row number $row."
    global npolymers nrows ncolumns
    if {$pID > $npolymers} {return}
    #if this is the first row, the column can be placed based on the first
    #polymer in the stack
    if {$row == 0} {
        puts stderr "Building to row starting from first stack polymer."
    } elseif { [expr $stack % 2] == 0 } {
        set pID [lay_brick $pID "+y"]
    } else {
        set pID [lay_brick $pID "\-y"]
    }

    for {set column 1} {$column < $ncolumns} {incr column} {
        if { [expr $row % 2] == [expr ($nrows*$stack)%2] } {
            set pID [lay_brick $pID "+x"]
        } else {
            set pID [lay_brick $pID "\-x" ]
        }
    }
    return $pID
}

proc lay_brick { padd direction } {
    global npolymers
    if {$padd > $npolymers} {return}
    puts stderr "Laying polymer brick number $padd in direction $direction."
    global brick_increment
    set lastP [atomselect top "segid P[expr $padd - 1]"]
    set poly [atomselect top "segid P$padd"]
    $poly moveby [vecinvert [measure center $poly]]
    $poly moveby [measure center $lastP]
    push_till_fit $padd [expr $padd - 1] $direction
    wiggle_till_fit $padd
    $lastP delete; $poly delete 
    return [expr $padd + 1]
}

proc push_till_fit { p1 p2 direction } {
    global brick_increment
    puts stderr "And puuuush!"
    set p1sel [atomselect top "segid P$p1"]
    set axis [string index $direction 1]
    set move_vec [get_axis_vec $axis]
    if {$axis == "z"} {
        set move_vec [vecscale 3.0 $move_vec]
    }
    set sign [string index $direction 0]
    if {$sign == "-"} {
        set move_vec [vecinvert $move_vec]
    }
    while {[brick_test $p1 $p2]} {
        $p1sel moveby [vecscale $brick_increment $move_vec]
    }
    $p1sel moveby [vecscale $brick_increment $move_vec]
    $p1sel delete
}

proc wiggle_till_fit { p } {
    global brick_increment
    puts stderr "Wiggle wiggle..."
    set psel [atomselect top "segid P$p"]
    while {[brick_test $p ""]} {
        $psel moveby [list \
            [expr $brick_increment*(2*rand()-1)] \
            [expr $brick_increment*(2*rand()-1)] \
            [expr $brick_increment*(0.5*rand()-0.25)] \
        ]
    }
    $psel delete
}

proc brick_test { p1 p2 } {
    global brick_increment
    if {$p2 == ""} {
        set testsel [atomselect top "(within $brick_increment of segid P$p1) and (not segid P$p1)"]
    } else {
        set testsel [atomselect top "segid P$p1 and within $brick_increment of segid P$p2"]
    }
    set testnum [$testsel num]
    $testsel delete;
    return $testnum
}

proc get_axis_vec { axis } {
    if {$axis == "x"} {
        return {1.0 0.0 0.0}
    }
    if {$axis == "y"} {
        return {0.0 1.0 0.0}
    } 
    if {$axis == "z"} {
        return {0.0 0.0 1.0}
    }
}

proc angry_baby { file } {
    global maximum_layer_throws maximum_polymer_throws
    #Throws the polymers in the file around until they don't clash
    mol new ${file}.psf waitfor all
    mol addfile ${file}.pdb
    set npolymers [llength [lsort -u [[atomselect top "all"] get segid]]]
    puts stderr "Throwing polymers around until they don't clash."
    for {set t 1} {$t <= [expr $maximum_layer_throws ]} {incr t} {
        puts -nonewline stderr "Layer throw $t  "
        throw_layer $npolymers
        set clashnum [count_all_clashes ]
        if { $clashnum == 0 } {
            puts stderr "Success. No clash." 
            break
        } else {
            puts stderr "has $clashnum clashing segments. Rethrowing single polymers."
        }
        puts stderr "Polymer ... fixed after ... throws"
        for {set p 1} {$p<=$npolymers} {incr p} {
            set pclash [count_segid_clashes P$p ]
            set s 0
            for {} {$s <= $maximum_polymer_throws && $pclash > 0} {incr s} {
                throw_polymer $npolymers P$p
                set pclash [count_segid_clashes P$p ]
            }
            if { $pclash > 0 } { 
                puts stderr "\nCouldn't fix clashes in segid P$p within $maximum_polymer_throws polymer throws."; 
                break 
            } else {
                puts -nonewline stderr "P$p:$s.." 
            }
        }
        if { $pclash == 0 } { 
            if { [count_all_clashes ] == 0 } { puts stderr "Success. No clash."; break }
        }
        if { $t == $maximum_layer_throws } {
            puts stderr "No success after $maximum_layer_throws tries. Creating new polymers." 
            return [build_layer $npolymers]
        }
    }
    puts stderr "\nWriting layer."

    set all [atomselect top "all"]
    $all moveby [vecinvert [measure center $all]]
    $all writepsf "angry_baby_build.psf"
    $all writepdb "angry_baby_build.pdb"
    return "angry_baby_build"
}

proc merge_polymers { polyfiles } {
    #Merge the polymers into a single pdb file
    package require topotools
    set midlist {}
    foreach poly $polyfiles {
        set mol [mol new $poly.psf waitfor all]
        mol addfile $poly.pdb
        lappend midlist $mol
    }
    set mol [::TopoTools::mergemols $midlist]
    animate write psf merged_polys.psf $mol
    animate write pdb merged_polys.pdb $mol
    mol delete all
    puts stderr "Polymers merged."
    return "merged_polys"
}

proc throw_layer { npolymers } {
    global space_per_polymer
    set box_size [expr $npolymers*$space_per_polymer]
    for {set p 1} {$p <= $npolymers} {incr p} {
        throw_polymer $box_size P$p
    }
}

proc throw_polymer { box_size segid } {
    set sel [atomselect top "segid $segid"]
    $sel moveby [vecinvert [measure center $sel]]
    set movebylist {}
    foreach axis {x y z} {
        #Random rotations
        $sel move [transaxis $axis [expr rand()*360]]
        #Random placement
        lappend movebylist [expr $box_size*rand()/2-$box_size/2]
    }
    $sel moveby $movebylist
    $sel delete
}

proc count_segid_clashes { segid } {
    global external_clash_tolerance
    set sel [atomselect top "(within $external_clash_tolerance of segid $segid) and (not segid $segid)"]
    set clashes [$sel num]
    $sel delete
    return $clashes
}

proc count_all_clashes {} {
    #Number of segids with clashes.
    global external_clash_tolerance
    set segidlist [lsort -u [[atomselect top "all"] get segid]]
    set clashsum 0
    foreach segid $segidlist {
        set sel [atomselect top "(within $external_clash_tolerance of segid $segid) and (not segid $segid)"]
        incr clashsum [llength [lsort -u [$sel get segid]]]      
    }
    return [expr $clashsum ]
}
