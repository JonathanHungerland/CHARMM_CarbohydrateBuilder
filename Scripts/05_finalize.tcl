proc finalize { file } {
    global simprepname
    set outfile [move_to_center $file]
    write_pbc $outfile
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

proc shift_stacks { file } {
   global nrows ncolumns nstacks
   mol new $file.psf
   mol addfile $file.pdb
   set all [atomselect top "all"]
   set mm [measure minmax $all]
   set x_length [expr abs([lindex $mm 0 0]-[lindex $mm 1 0])]
   set y_length [expr abs([lindex $mm 0 1]-[lindex $mm 1 1])]
   set polys_per_stack [expr $nrows*$ncolumns]
   for {set s 0} {$s <= $nstacks} {incr s} {
      set seltext "segid "
      for {set p 1} {$p<=$polys_per_stack} {incr p} {
         set add_poly [expr $s*$polys_per_stack+$p]
         set seltext "${seltext} P${add_poly}"
      }
      set sel [atomselect top $seltext]
      $sel moveby [list \
                   [expr $x_length*(rand()-0.5)] \
                   [expr $y_length*(rand()-0.5)] \
                   0]
   }
   $all writepdb shifted.pdb
   $all writepsf shifted.psf
   return "shifted"
}
