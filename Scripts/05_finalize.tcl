proc finalize { file } {
    global simprepname
    set outfile [move_to_center $file]
    if { $build_method == "bored_bricklayer" } {
        
    } else { write_pbc $outfile }
    return $outfile
}

proc move_to_center { file } {
    mol new $file.psf
    mol addfile $file.pdb
    set all [atomselect top "all"]
    $all moveby [vecinvert [measure center $all]]
    $all writepsf centered.psf
    $all writepdb centered.pdb
    $all delete
    return "centered"
}

proc write_pbc { file } {
    mol new $file.psf
    mol addfile $file.pdb
    set all [atomselect top "all"]
    #pbc output
    $all moveby [vecinvert [measure center $all]]
    set mm [measure minmax $all]
    set vec [vecsub [lindex $mm 1] [lindex $mm 0]]
    set pbcfile [open "${file}.pbc" w]
    puts $pbcfile "xpbc=[lindex $vec 0]"
    puts $pbcfile "ypbc=[lindex $vec 1]"
    puts $pbcfile "zpbc=[lindex $vec 2]"
    close $pbcfile
}

