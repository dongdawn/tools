# 20111223 V2.2: faster by defining more and more detailed selections #
# 20111223 V2.1: faster by defining more detailed selections #
#                sel1 and sel2 must IDENTICAL or ABSULUTELY DIFFERENT, or there will be error caused by VMD
# 20111221 V2.0: seltxt1 & seltxt2 can be the same #
proc HBpairAnal {outfile nHIS HISframe HBprcnt molid HBlength HBangle seltxt1 seltxt2} {
### option ###
### 所有帧的分段个数nHIS，出现频率底限HBprcnt ###
### 其中，hbABnum为A-H-B氢键总数，hbListBA为B-H-A氢键总数
###				HBList($unt)为$unt原子涉及氢键的次数
###				HBFrmList($unt$cF$iHIS)为第iHIS时间段内$unt原子涉及氢键的次数
###				HBpairList($hbH$cA$hbA)为某氢键出现的总次数
###				HBpairFrmList($hbH$cA$hbA$cF$iHIS)为第iHIS时间段内某氢键出现的次数
## parameters ##
  set outFL [open $outfile w]
  set outHISFL [open ${outfile}HIS w]
  set outHISprFL [open ${outfile}HISpair w]
  set outHISprPFL [open ${outfile}HISpairP w]
  set nframe [molinfo $molid get numframes]
 #  set HISframe [expr "$nframe/$nHIS"]
  set HBlow [expr "int($nframe*$HBprcnt)"]
  set selALL [atomselect $molid all]
  set cA A
  set cF F
  if {[info exists HBList]} {unset HBList}
  if {[info exists HBFrmList]} {unset HBFrmList}
  if {[info exists HBpairList]} {unset HBpairList}
  if {[info exists HBpairFrmList]} {unset HBpairFrmList}
## set selList ##
  if {$seltxt1 != $seltxt2} then {
  	set resSEList [atomselect $molid "$seltxt1 and name \"O.*\" \"N.*\""]
  	set dmnSEList [atomselect $molid "$seltxt2 and name \"O.*\" \"N.*\""]
  } else {
    set resSEList ""
    set dmnSEList ""
    set dtilSeList {"name \"O.*\""  "name \"N.*\""                 \
                    "name O"        "name \"O.*\" and not name O"  \
                    "name N"        "name \"N.*\" and not name N"  \
                   }
    foreach {ur ud} $dtilSeList {
      set selA [atomselect $molid "$seltxt1 and $ur"]
      set selB [atomselect $molid "$seltxt2 and $ud"]
      if {[$selA num]>0 && [$selB num]>0} then {
        lappend resSEList $selA
        lappend dmnSEList $selB
      }
    }
  # 
    set snglSeList {"OD1 OD2 OE1 OE2 OG OG1 OH OT1 OT2"   "ND1 ND2 NE NE1 NE2 NH1 NH2 NZ"}
    foreach uss $snglSeList {
      set ndx_1 [expr [llength $uss] -1]
      for {set idx 0} {$idx < $ndx_1} {incr idx} {
        set selA [atomselect $molid "$seltxt1 and name [lindex $uss $idx]"]
        set selB [atomselect $molid "$seltxt1 and name [lrange $uss [expr $idx +1] $ndx_1]"]
        if {[$selA num]>0 && [$selB num]>0} then {
          lappend resSEList $selA
          lappend dmnSEList $selB
        } 
      }
    #  
      foreach us $uss {
        set seltp [atomselect $molid "$seltxt1 and name $us"]
        if {[$seltp num]>1} then {
        set idces [$seltp get index]
        set ndx_1 [expr [$seltp num] -1]
        for {set idx 0} {$idx < $ndx_1} {incr idx} {
          lappend resSEList [atomselect $molid "index [lindex $idces $idx]"]
          lappend dmnSEList [atomselect $molid "index [lrange $idces [expr $idx +1] $ndx_1]"]
        }
        }
      }
    }
  }
## loop in frames ##
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
 # overall summany #
    set hbsum 0
    set hbABsum 0
    set hbBAsum 0
    set iHIS [expr "int($iframe/$HISframe)"]
    foreach urs $resSEList uds $dmnSEList {
      $urs frame $iframe
      $uds frame $iframe
      set hbListAB [measure hbonds $HBlength $HBangle $urs $uds]
      set hbListBA [measure hbonds $HBlength $HBangle $uds $urs]
      set hbABnum [llength [lindex $hbListAB 0]]
      set hbBAnum [llength [lindex $hbListBA 0]]
      set hbABsum [expr $hbABsum + $hbABnum]
      set hbBAsum [expr $hbBAsum + $hbBAnum]
      set hbsum [expr $hbsum + $hbABnum + $hbBAnum]
 # atom details #
      foreach hbA [lindex $hbListAB 1] hbH [lindex $hbListAB 2] {
        if {![info exists HBList($hbA)]} {set HBList($hbA) 0}
        incr HBList($hbA)
        if {![info exists HBList($hbH)]} {set HBList($hbH) 0}
        incr HBList($hbH)
        if {![info exists HBFrmList($hbH$cF$iHIS)]} {set HBFrmList($hbH$cF$iHIS) 0}
        incr HBFrmList($hbH$cF$iHIS)
        if {![info exists HBFrmList($hbA$cF$iHIS)]} {set HBFrmList($hbA$cF$iHIS) 0}
        incr HBFrmList($hbA$cF$iHIS)
        if {![info exists HBpairList($hbH$cA$hbA)]} {set HBpairList($hbH$cA$hbA) 0}
        incr HBpairList($hbH$cA$hbA)
        if {![info exists HBpairFrmList($hbH$cA$hbA$cF$iHIS)]} {set HBpairFrmList($hbH$cA$hbA$cF$iHIS) 0}
        incr HBpairFrmList($hbH$cA$hbA$cF$iHIS)
      }
      foreach hbA [lindex $hbListBA 1] hbH [lindex $hbListBA 2] {
        if {![info exists HBList($hbA)]} {set HBList($hbA) 0}
        incr HBList($hbA)
        if {![info exists HBList($hbH)]} {set HBList($hbH) 0}
        incr HBList($hbH)
        if {![info exists HBFrmList($hbA$cF$iHIS)]} {set HBFrmList($hbA$cF$iHIS) 0}
        incr HBFrmList($hbA$cF$iHIS)
        if {![info exists HBFrmList($hbH$cF$iHIS)]} {set HBFrmList($hbH$cF$iHIS) 0}
        incr HBFrmList($hbH$cF$iHIS)
        if {![info exists HBpairList($hbA$cA$hbH)]} {set HBpairList($hbA$cA$hbH) 0}
        incr HBpairList($hbA$cA$hbH)
        if {![info exists HBpairFrmList($hbA$cA$hbH$cF$iHIS)]} {set HBpairFrmList($hbA$cA$hbH$cF$iHIS) 0}
        incr HBpairFrmList($hbA$cA$hbH$cF$iHIS)
      }
    }
    puts $outFL "$iframe $hbABsum $hbBAsum $hbsum"
  }
  
## puts HBList ##
  foreach {nm vl} [array get HBList] {
    set outtxt ""
    if {$vl > $HBlow} {
      lappend outtxt $nm
      lappend outtxt [lindex [$selALL get chain] $nm]
      lappend outtxt [lindex [$selALL get resname] $nm]
      lappend outtxt [lindex [$selALL get resid] $nm]
      lappend outtxt [lindex [$selALL get name] $nm]
      lappend outtxt $vl
      lappend outtxt :
      for {set iHIS 0} {$iHIS < $nHIS} {incr iHIS} {
        if {![info exists HBFrmList($nm$cF$iHIS)]} {set HBFrmList($nm$cF$iHIS) 0}
        lappend outtxt $HBFrmList($nm$cF$iHIS)
      }
      puts $outHISFL $outtxt
    }
  }
## puts HBpairList ##
  foreach {nm vl} [array get HBpairList] {
    set outtxt ""
    set whA [string first $cA $nm]
    set nm1 [string range $nm 0 [expr $whA-1]]
    set nm2 [string range $nm [expr $whA+1] [expr [string length $nm]-1]]
    if {$vl > $HBlow} {
      lappend outtxt $nm1
      lappend outtxt [lindex [$selALL get chain] $nm1]
      lappend outtxt [lindex [$selALL get resname] $nm1]
      lappend outtxt [lindex [$selALL get resid] $nm1]
      lappend outtxt [lindex [$selALL get name] $nm1]
      lappend outtxt -
      lappend outtxt $nm2
      lappend outtxt [lindex [$selALL get chain] $nm2]
      lappend outtxt [lindex [$selALL get resname] $nm2]
      lappend outtxt [lindex [$selALL get resid] $nm2]
      lappend outtxt [lindex [$selALL get name] $nm2]
      set outtxtp $outtxt 
      puts $outHISprPFL "$outtxtp [expr $vl / double($nframe)]"
      lappend outtxt $vl
      lappend outtxt :
      for {set iHIS 0} {$iHIS < $nHIS} {incr iHIS} {
        if {![info exists HBpairFrmList($nm1$cA$nm2$cF$iHIS)]} {set HBpairFrmList($nm1$cA$nm2$cF$iHIS) 0}
        lappend outtxt $HBpairFrmList($nm1$cA$nm2$cF$iHIS)
      }
      puts $outHISprFL $outtxt
    }
  }
  close $outFL
  close $outHISFL
  close $outHISprFL
  close $outHISprPFL
}
# outfile nHIS HISframe HBprcnt molid HBlength HBangle seltxt1 seltxt2

