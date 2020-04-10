proc Smallmolxyz {outfile molid sel1} {
  set fl [open $outfile w]
  set motif1 [atomselect $molid "$sel1"]
  set nframe [molinfo $molid get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $motif1 frame $iframe
    puts $fl "$iframe [measure center $motif1 weight mass]"
  }
  close $fl
}

#Smallmolxyz resid506 12 "resid 506"

#statistics small molecular in certain x y of protein 
proc Smallmolxy_in {outfile molid sel1} {
  set fl [open $outfile w]
  set motif1 [atomselect $molid "$sel1"]
  set mima [measure minmax [atomselect $molid "protein"]]
  set minx [lindex [lindex $mima 0] 0]
  set maxx [lindex [lindex $mima 1] 0]
  set miny [lindex [lindex $mima 0] 1]
  set maxy [lindex [lindex $mima 1] 1]
  set nframe [molinfo $molid get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $motif1 frame $iframe
    set molxyz [measure center $motif1 weight mass]
    if {[lindex $molxyz 0] >$minx && [lindex $molxyz 0] <$maxx && [lindex $molxyz 1] > $miny && [lindex $molxyz 1] < $maxy} then {
    	
    	puts $fl "$iframe [lindex $molxyz 2]"
    	}
    
  }
  close $fl
}

#Smallmolxy_in resid500 12 "resid 500"
for {set i 500} {$i < 512} {incr i} {
	Smallmolxy_in resid${i} 0 "resid $i"
}


#statistics small mol in each channel   need to define z by yourself  

proc molinChan {outfile molid sel1} {
  set fl [open $outfile w]
  set motif1 [atomselect $molid "$sel1"]

  set minz 29
  set maxz 62
  set nframe [molinfo $molid get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $motif1 frame $iframe
    set mima [measure minmax [atomselect $molid "protein"]]
    set minx [lindex [lindex $mima 0] 0]
    set maxx [lindex [lindex $mima 1] 0]
    set miny [lindex [lindex $mima 0] 1]
    set maxy [lindex [lindex $mima 1] 1]
    set molxyz [measure center $motif1 weight mass]
    if {[lindex $molxyz 0] > $minx && [lindex $molxyz 0] < 53 && [lindex $molxyz 1] > $miny && [lindex $molxyz 1] < 53 && [lindex $molxyz 2] > $minz && [lindex $molxyz 2] <$maxz} {
    	set pore 1
    } elseif {[lindex $molxyz 0] > $minx && [lindex $molxyz 0] < 53 && [lindex $molxyz 1] < $maxy && [lindex $molxyz 1] > 53 && [lindex $molxyz 2] > $minz && [lindex $molxyz 2] <$maxz} {
    	set pore 2
    } elseif {[lindex $molxyz 0] < $maxx && [lindex $molxyz 0] > 53 && [lindex $molxyz 1] < $maxy && [lindex $molxyz 1] > 53 && [lindex $molxyz 2] > $minz && [lindex $molxyz 2] <$maxz} {
    	set pore 3
    } elseif {[lindex $molxyz 0] < $maxx && [lindex $molxyz 0] > 53 && [lindex $molxyz 1] > $miny && [lindex $molxyz 1] < 53 && [lindex $molxyz 2] > $minz && [lindex $molxyz 2] <$maxz} {
    	set pore 4
    } else {
    	set pore 0
    }
    puts $fl "$iframe $pore"
  }
  close $fl
}


#molinChan resid500inpore 12 "resid 500"

#for {set i 500} {$i < 512} {incr i} {
#	molinChan resid${i}inpore  12 "resid $i"
#}

