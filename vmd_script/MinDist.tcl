proc min {list} {
  set min 100.0
  foreach unt $list {
    if {$unt < $min} then {set min $unt}
  }
  return $min
}
proc MinDist {outfile molid1 sel1 molid2 sel2} {
  set fl [open $outfile w]
  set motif1 [atomselect $molid1 "$sel1"]
  set motif2 [atomselect $molid2 "$sel2"]
  set nframe [molinfo $molid1 get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $motif1 frame $iframe
    $motif2 frame $iframe
    set distList ""
    foreach xyz1 [$motif1 get {x y z}] {
      foreach xyz2 [$motif2 get {x y z}] {
        lappend distList [vecdist $xyz1 $xyz2]
      }
    }
  puts $fl "$iframe [min $distList]"
  }
  close $fl
}

MinDist mdist_301312c 0 "chain C and resid 301 and sidechain and noh" 0 "chain C and resid 312 and sidechain and noh"
MinDist mdist_301312d 0 "chain D and resid 301 and sidechain and noh" 0 "chain D and resid 312 and sidechain and noh"




# MalFGK2-MBP
MinDist dist_MalK_A165-B38 0 "chain A and resid 165" 0 "chain B and resid 38"
MinDist dist_MalK_B165-A38 0 "chain B and resid 165" 0 "chain A and resid 38"
MinDist dist_MalK_A166-B192 0 "chain A and resid 166" 0 "chain B and resid 192"
MinDist dist_MalK_B166-A192 0 "chain B and resid 166" 0 "chain A and resid 192"



MinDist dist_pins-insc 0 "chain A and resid 150 to 160" 3 "chain B"

# BtuCD
MinDist dist_AB143m 0 "chain A and resid 143" "chain B and resid 143"
MinDist dist_AB146m 0 "chain A and resid 146" "chain B and resid 146"
MinDist dist_AB147m 0 "chain A and resid 147" "chain B and resid 147"
MinDist dist_AB90m 0 "chain A and resid 90" "chain B and resid 90"
MinDist dist_AB82-89m 0 "chain A and resid 82 to 89" "chain B and resid 82 to 89"


