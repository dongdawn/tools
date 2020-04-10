#####################
# draw_bond function#
#####################

#Usage:
# draw_bond -mol1 0 -index1 1 -mol2 1 -index2 1 -h_nbars 25 -h_space 1 -h_arrow 0 -h_radius 0.1 -color yellow -mat Opaque -h_type cylinder -h_resol 25
proc draw_bond_help { } {
	puts "----------draw bond function help information-------------
	#Usage example:
	#draw_bond -mol1 0 -index1 1 -mol2 1 -index2 1 -h_nbars 25 -h_space 1 -h_arrow 0 -h_radius 0.1 -color yellow -mat Opaque -h_type cylinder -h_resol 25
	#--->options:
	# -mol1 0; the first molecule id, the default value is top
	# -index1 200; the first atom index number in mol1
	# -mol2 0; the the molecular id of the second atom
	# -index2; the second atom index number in mol2
	# -h_nbars 25; number of bars
	# -h_space 1; space between bars
	# -h_arrow 0; draw arrows on the two sides of the bond or not
	# -h_radius 0.1; the bond radius
	# -color yellow; the bond color
	# -mat Opaque; the material render method
	# -h_type cylinder; supporting sphere, cylinder, pymol, cone, line
	# -h_reol 25; bond resolution"
}

proc draw_bond {args} {
	#the original function is: draw_bond {res1 res2 {color yellow}}
	#set default values
	#set mol_graphics_id 0

	set mol_id_1 top
	set mol_id_2 top

	set hbonds_nbars 25
	set hbonds_ratio 1
	set hbonds_arrow 0
	set hbonds_radius 0.1


	set color yellow
	set material Opaque
	set h_resol 25

	set h_type cylinder

	#Parse options---------------------------
	# draw_bond -mol1 0 -index1 1 -mol2 1 -index2 1 -h_nbars 25 -h_space 1 -h_arrow 0 -h_radius 0.1 -color yellow -h_type cylinder
	for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {
		set arg [ lindex $args $argnum ]
		set val [ lindex $args [expr $argnum + 1]]

		switch -- $arg {
	    	"-mol1"      {set mol_id_1 $val; incr argnum}
	    	"-index1"	 {set index1 $val; incr argnum}
	    	"-mol2"	     {set mol_id_2 $val; incr argnum}
	    	"-index2"	 {set index2 $val; incr argnum}
	    	"-h_nbars"   {set hbonds_nbars $val; incr argnum}
	    	"-h_space"   {set hbonds_ratio $val; incr argnum}
	    	"-h_radius"  {set hbonds_radius $val; incr argnum}
	    	"-h_arrow"   {set hbonds_arrow $val; incr argnum}
	    	"-color"     {set color $val; incr argnum}
	    	"-mat"       {set material $val; incr argnum}
	    	"-h_resol"   {set h_resol $val; incr argnum}
	    	"-h_type"    {set h_type $val; incr argnum}

	    	default { puts "unknown option: $arg and you must set the -index1 and index2 values "; return }
		}
    }
	#-----------stop parse-------------------
	lassign [[atomselect $mol_id_1 "index $index1"] get {x y z}] START
	lassign [[atomselect $mol_id_2 "index $index2"] get {x y z}] END
	set nbars $hbonds_nbars
	set X1 [lindex $START 0] ; set X2 [lindex $END 0]
	set Y1 [lindex $START 1] ; set Y2 [lindex $END 1]
	set Z1 [lindex $START 2] ; set Z2 [lindex $END 2]
	if { $X1 != "" && $X2 != "" && $Y1 != "" && $Y2 != "" && $Z1 != "" && $Z2 != "" } {
		set drawarrows $hbonds_arrow
		set LX [expr ( $X2 - $X1 )]
		set LY [expr ( $Y2 - $Y1 )]
		set LZ [expr ( $Z2 - $Z1 )]

		if { $drawarrows != 0 } {
			incr nbars 2
			set imin 1
		} else {
			set imin 0
		}




		#graphics $mol_graphics_id color $color
		draw color $color
		draw material $material
		set line_radius [expr ($hbonds_radius * 20)]


		set barx [expr ($hbonds_ratio * $LX) / ($nbars * ($hbonds_ratio + 1) - 1)]
		set bary [expr ($hbonds_ratio * $LY) / ($nbars * ($hbonds_ratio + 1) - 1)]
		set barz [expr ($hbonds_ratio * $LZ) / ($nbars * ($hbonds_ratio + 1) - 1)]

		set gapx [expr ($barx / $hbonds_ratio)]
		set gapy [expr ($bary / $hbonds_ratio)]
		set gapz [expr ($barz / $hbonds_ratio)]

		for {set i $imin} { $i < [expr $nbars - $imin]} { incr i 1} {
			set DEBX [expr ($X1 + ($i * ($barx + $gapx)))]
			set ENDX [expr ($DEBX + $barx)]
			set DEBY [expr ($Y1 + ($i * ($bary + $gapy)))]
			set ENDY [expr ($DEBY + $bary)]
			set DEBZ [expr ($Z1 + ($i * ($barz + $gapz)))]
			set ENDZ [expr ($DEBZ + $barz)]

			if { $h_type == "sphere" } {
				draw sphere "$DEBX $DEBY $DEBZ" resolution $h_resol radius $hbonds_radius
			} elseif {$h_type == "cone"} {
				draw cone "$DEBX $DEBY $DEBZ" "$ENDX $ENDY $ENDZ" resolution $h_resol radius $hbonds_radius
			} elseif {$h_type == "line"} {
				draw line $START $END style dashed width $line_radius
			} elseif {$h_type == "pymol"} {
				draw cylinder "$DEBX $DEBY $DEBZ" "$ENDX $ENDY $ENDZ" resolution $h_resol radius $hbonds_radius filled yes
				draw sphere "$DEBX $DEBY $DEBZ" resolution $h_resol radius $hbonds_radius
				draw sphere "$ENDX $ENDY $ENDZ" resolution $h_resol radius $hbonds_radius
			} else {
				draw cylinder "$DEBX $DEBY $DEBZ" "$ENDX $ENDY $ENDZ" resolution $h_resol radius $hbonds_radius filled yes
			}


			if { $i eq 1 && $drawarrows } {
				lassign "$DEBX $DEBY $DEBZ" ARROWHEADX ARROWHEADY ARROWHEADZ
			}
			if { $i eq [expr $nbars - 2] && $drawarrows } {
				lassign "$ENDX $ENDY $ENDZ" ARROWTAILX ARROWTAILY ARROWTAILZ
			}
			unset DEBX DEBY DEBZ ENDX ENDY ENDZ
		}

		if { $drawarrows == 2} {

			draw cone "$ARROWHEADX $ARROWHEADY $ARROWHEADZ" "$X1 $Y1 $Z1" radius $hbonds_radius resolution $h_resol
			draw cone "$ARROWTAILX $ARROWTAILY $ARROWTAILZ" "$X2 $Y2 $Z2" radius $hbonds_radius resolution $h_resol

			unset ARROWHEADX ARROWHEADY ARROWHEADZ ARROWTAILX ARROWTAILY ARROWTAILZ
		} elseif { $drawarrows == 1} {
            draw cone "$ARROWHEADX $ARROWHEADY $ARROWHEADZ" "$X1 $Y1 $Z1" radius $hbonds_radius resolution $h_resol
            unset ARROWHEADX ARROWHEADY ARROWHEADZ ARROWTAILX ARROWTAILY ARROWTAILZ
		} elseif { $drawarrows == -1} {
            draw cone "$ARROWTAILX $ARROWTAILY $ARROWTAILZ" "$X2 $Y2 $Z2" radius $hbonds_radius resolution $h_resol
			unset ARROWHEADX ARROWHEADY ARROWHEADZ ARROWTAILX ARROWTAILY ARROWTAILZ
		}
		unset X1 X2 Y1 Y2 Z1 Z2
		unset LX LY LZ barx bary barz gapx gapy gapz line_radius
	}
}

draw_bond -mol1 0 -index1 3009 -mol2 0 -index2 3870 -h_nbars 15 -h_space 1 -h_arrow 0 -h_radius 0.1 -color silver -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 3003 -mol2 0 -index2 3870 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.1 -color gray -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 2903 -mol2 0 -index2 3871 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.1 -color gray -mat Opaque -h_type cylinder -h_resol 10


draw_bond -mol1 0 -index1 3009 -mol2 0 -index2 3949 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.1 -color gray -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 3003 -mol2 0 -index2 3949 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.1 -color silver -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 2903 -mol2 0 -index2 3950 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.1 -color black -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 2923 -mol2 0 -index2 3945 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.1 -color black -mat Opaque -h_type cylinder -h_resol 10

draw_bond -mol1 0 -index1 2903 -mol2 0 -index2 3945 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.09 -color black -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 2923 -mol2 0 -index2 3941 -h_nbars 10 -h_space 1.6 -h_arrow 0 -h_radius 0.08 -color yellow -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 3003 -mol2 0 -index2 3944 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.09 -color black -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 3009 -mol2 0 -index2 3949 -h_nbars 10 -h_space 1.2 -h_arrow 0 -h_radius 0.08 -color pink -mat Opaque -h_type cylinder -h_resol 10

draw_bond -mol1 0 -index1 2903 -mol2 0 -index2 3973 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.1 -color green -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 2903 -mol2 0 -index2 3969 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.1 -color orange -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 2923 -mol2 0 -index2 3941 -h_nbars 10 -h_space 1.6 -h_arrow 0 -h_radius 0.09 -color gray -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 3003 -mol2 0 -index2 3944 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.1 -color black -mat Opaque -h_type cylinder -h_resol 10
draw_bond -mol1 0 -index1 3009 -mol2 0 -index2 3968 -h_nbars 10 -h_space 1 -h_arrow 0 -h_radius 0.09 -color yellow -mat Opaque -h_type cylinder -h_resol 10

draw_bond -mol1 32 -index1 210 -mol2 32 -index2 2820 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25
draw_bond -mol1 32 -index1 469 -mol2 32 -index2 2820 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25
draw_bond -mol1 32 -index1 837 -mol2 32 -index2 2820 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25
draw_bond -mol1 32 -index1 1002 -mol2 32 -index2 2820 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25
draw_bond -mol1 32 -index1 210 -mol2 32 -index2 3034 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25
draw_bond -mol1 32 -index1 469 -mol2 32 -index2 3034 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25
draw_bond -mol1 32 -index1 837 -mol2 32 -index2 3034 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25
draw_bond -mol1 32 -index1 1002 -mol2 32 -index2 3034 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25
draw_bond -mol1 32 -index1 2484 -mol2 32 -index2 2648 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25
draw_bond -mol1 32 -index1 2484 -mol2 32 -index2 3034 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25
draw_bond -mol1 32 -index1 2648 -mol2 32 -index2 3034 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color red -mat Opaque -h_type cylinder -h_resol 25

draw_bond -mol1 32 -index1 1002 -mol2 32 -index2 3034 -h_nbars 25 -h_space 2 -h_arrow 0 -h_radius 0.2 -color black -mat Opaque -h_type cylinder -h_resol 25