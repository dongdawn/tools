proc Ax-Ax {outfile molID selax1txt selax2txt} {
  set fl [open $outfile w]
  set selax1 [atomselect $molID "$selax1txt"]
  set selax2 [atomselect $molID "$selax2txt"]
  set nframe [molinfo $molID get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $selax1 frame $iframe
    $selax2 frame $iframe
    set ax1 [lindex [lindex [measure inertia $selax1] 1] 2]
    set ax2 [lindex [lindex [measure inertia $selax2] 1] 2]
    set tp [expr acos([vecdot $ax1 $ax2])/3.1416*180.0]
    puts $fl "$iframe $tp"
  }
  close $fl
}
#Ax-Ax ax-ax 0 "chain A and resid 142 to 167 and name CA" "chain B and resid 142 to 167 and name CA"