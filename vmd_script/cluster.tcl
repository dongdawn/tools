### calculate cluster dcds ###
# cutoff #
foreach u {0.01 0.1 0.5} {
  set wfl [open clu${u}.log w]
  puts $wfl [measure cluster [atomselect 8 "resname GOL and noh"] distfunc rmsd weight mass cutoff $u]
  close $wfl
}

foreach u {1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5} {
  set wfl [open clu${u}.log w]
  puts $wfl [measure cluster [atomselect 0 "resname FAR and noh"] num 20 distfunc rmsd weight mass cutoff $u]
  close $wfl
}

foreach u {3.0} {
  set wfl [open clu${u}.log w]
  puts $wfl [measure cluster [atomselect 0 "protein and backbone"] num 100 distfunc rmsd weight mass cutoff $u]
  close $wfl
}

### cluster summary ###
#foreach u {1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0} {
foreach u {3.0} {
  set rfl [open clu${u}.log r]
  set wfl [open clu${u}.txt w]
  gets $rfl line
  close $rfl
  puts "### clu${u}.log ###"
  set nclu [llength $line]
  for {set iclu 0} {$iclu < $nclu} {incr iclu} {
  	puts $wfl [lindex $line $iclu]
    #puts "# [expr $iclu +1] # middle=[lindex [lindex $line $iclu] 0] num=[llength [lindex $line $iclu]]"
  }
  close $wfl
  close $rfl
}


### produce cluster timeline ###
set nclu 4
set logfile clu2.0.log
set tlfile clusterline 
# clu2.0.log for far
# clu3.2.log for dpc

set rfl [open $logfile r]
gets $rfl line
close $rfl

set wfl [open $tlfile w]
set oclu [expr $nclu + 1]

set iclu 0
foreach clulist $line {
  incr iclu
  if {$iclu <= $nclu} then {
    foreach u $clulist {set Cndx($u) $iclu}
  } else {
    foreach u $clulist {set Cndx($u) $oclu}
  }
}

foreach u [lsort -integer [array names Cndx]] {
  puts $wfl "$u $Cndx($u)"
}
close $wfl




### produce cluster dcds ###
set nclu 5
set pdbfile ../prot.pdb
set dcdfile 9k.dcd
set logfile clu3.2.log
# clu2.0.log for far
# clu3.2.log for dpc

set rfl [open $logfile r]
gets $rfl line
close $rfl

for {set iclu 1} {$iclu <= $nclu} {incr iclu} {
  set center [lindex [lindex $line [expr $iclu -1]] 0]
  set molid [mol new $pdbfile]
  mol addfile $dcdfile type dcd first $center last $center waitfor all
  animate delete beg 0 end 0
  animate write pdb clu${iclu}.pdb $molid
  mol delete $molid
}

for {set iclu 1} {$iclu <= $nclu} {incr iclu} {
  set molid [mol new $pdbfile]
  mol addfile $dcdfile type dcd waitfor all
  animate delete beg 0 end 0
  set frmlist [lsort -integer [lindex $line [expr $iclu -1]]]
  set ifrm 0
  set fstfrm [lindex $frmlist 0]
  foreach u $frmlist {
    incr ifrm
    set lastfrm [lindex $frmlist [expr [llength $frmlist] -1]]
    set nfrm [expr [molinfo $molid get numframes] - $ifrm]
    if {$lastfrm < $nfrm} then {
      animate delete beg [expr $lastfrm +1] end $nfrm $molid
    }
    set frmlist [lrange $frmlist 0 [expr [llength $frmlist] -2]]
  }
  if {$fstfrm != 0} then {
    animate delete beg 0 end [expr $fstfrm -1] $molid
  }
  animate write dcd clu${iclu}.dcd $molid
  mol delete $molid
}




