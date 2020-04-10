proc MCenterDist {outfile molid sel1 sel2} {
  set fl [open $outfile w]
  set motif1 [atomselect $molid "$sel1"]
  set motif2 [atomselect $molid "$sel2"]
  set nframe [molinfo $molid get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $motif1 frame $iframe
    $motif2 frame $iframe
    puts $fl "$iframe [vecdist [measure center $motif1 weight mass] [measure center $motif2 weight mass]]"
  }
  close $fl
}