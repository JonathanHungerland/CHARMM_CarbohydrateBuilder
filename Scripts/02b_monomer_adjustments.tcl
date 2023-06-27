proc define_plane { sel1 sel2 sel3 } {
  #defines a plane based on the center of three selections
  #returns then normal vector of the plane
  set p1 [measure center [atomselect top "$sel1"]]
  set p2 [measure center [atomselect top "$sel2"]]
  set p3 [measure center [atomselect top "$sel3"]]
  set a [vecsub $p2 $p1]
  set b [vecsub $p3 $p1]
  return [veccross $a $b]
}

proc measure_angle { v axis } {
  #returns the angle between two vectors
  set L [veclength $v]
  set dot [vecdot $v [axis_vector $axis]]
  return [expr acos($dot/$L)]
}

proc axis_vector { axis } {
  if { $axis == "x" } {
     return [list 1 0 0]
  }
  if { $axis == "y" } {
     return [list 0 1 0]
  }
  return [list 0 0 1]
}

proc adjust_plane_to_axis { move_sel plane_sel1 plane_sel2 plane_sel3 axis } {
  set normal [define_plane $plane_sel1 $plane_sel2 $plane_sel3]
  set m_sel [atomselect top "$move_sel"]
  $m_sel move [transvecinv $normal]
  foreach a [list "x" "y" "z"] {
    if { $a == $axis } {
      continue
    }
    $m_sel move [transaxis $a 90]
  }
  set angle [measure_angle $normal $axis]
  puts "After rotation angle of $angle degrees between normal and axis $axis"
}

proc center_sel { move_sel ref_sel } {
  set ref [atomselect top "$ref_sel"]
  set move [atomselect top "$move_sel"]
  $move moveby [vecinvert [measure center $ref]]
}
