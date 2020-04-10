proc frespair {sel1text sel2text cutoff outfile} {

	set outputs [open "$outfile" w]

	#########################################################
	if {[info exists atomIDarr1]} {unset atomIDarr1}
	if {[info exists atomIDarr2]} {unset atomIDarr2}
	if {[info exists prcnt]}      {unset prcnt}
	if {[info exists distarr]}    {unset distarr}

	set molid [molinfo top]
	set sel1  [atomselect $molid "$sel1text and name CA"]
	set sel2  [atomselect $molid "$sel2text and name CA"]
	set sel1resIDList [lsort -integer -increasing [$sel1 get resid]]
	set sel2resIDList [lsort -integer -decreasing [$sel2 get resid]]

	$sel1 delete
	$sel2 delete
	
	foreach sel1resID $sel1resIDList {
		set sel1res [atomselect $molid "resid $sel1resID"]
		set atomIDarr1($sel1resID) [$sel1res get index]
		$sel1res delete
	}
	foreach sel2resID $sel2resIDList {
		set sel2res [atomselect $molid "resid $sel2resID"]
		set atomIDarr2($sel2resID) [$sel2res get index]
		$sel2res delete
	}
  
	#########################################################
	set totF [molinfo $molid get numframes]
	
	foreach sel2resID $sel2resIDList {
		set lineOut ""
		append lineOut [format "%3i:  " $sel2resID]
		foreach sel1resID $sel1resIDList {
			
			if {[info exists prcnt($sel1resID:$sel2resID)]} {
				append lineOut [format "%5.1f " $prcnt($sel1resID:$sel2resID)]
			} else {
				for {set i 0} {$i<$totF} {incr i} {
					set distList ""
					foreach sel1AtomID $atomIDarr1($sel1resID) {
						foreach sel2AtomID $atomIDarr2($sel2resID) { 
							if {$sel1AtomID == $sel2AtomID} {
								lappend distList "0.0"
							} else {
								lappend distList [measure bond "$sel1AtomID $sel2AtomID" molid $molid frame $i]
							}
						}
					}
					set miniDist [lindex [lsort -real -increasing $distList] 0]
					if {$miniDist<=$cutoff} {
						set distarr($sel2resID:$sel1resID@$i) 1
					} else {
						set distarr($sel2resID:$sel1resID@$i) 0
					}
				}

				set count 0;
				for {set i 0} {$i<$totF} {incr i} {
					if {$distarr($sel2resID:$sel1resID@$i)>0} {
						incr count
					} 
				}
				set prcnt($sel2resID:$sel1resID) [expr $count*100/($totF+0.0)]
				append lineOut [format "%5.1f " $prcnt($sel2resID:$sel1resID)]
			}
		}

		puts $outputs $lineOut
#		flush $outputs
	}
	close $outputs

}

#############################################
#set sel1text "chain B and switchloop and noh";
#set sel2text "chain B and switchloop and noh";
set sel1text "(resid 1 to 67) and noh";
set sel2text "(resid 1 to 67) and noh";
set cutoff   "3.0";
set outfile  "out.dat"
#############################################

frespair $sel1text $sel2text $cutoff $outfile

#exit 0;


