#count the contact 
proc contactCountAtom {outfile molID cutoff selax1txt selax2txt} {
  set fl [open $outfile w]
  set sel1 [atomselect $molID "$selax1txt"]
  set sel2 [atomselect $molID "$selax2txt"]
  set nframe [molinfo $molID get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $sel1 frame $iframe
    $sel2 frame $iframe
    set contactNum [llength [lindex [measure contacts $cutoff $sel1 $sel2] 1]]
    puts $fl "$iframe        $contactNum"
  }
  close $fl
  puts "done"
}
#contactCountAtom contactNumber 0 3.5 "chain A " "chain B"
#llength [lindex [measure contacts 3.5 [atomselect top "chain A"] [atomselect top "chain B"]] 1]