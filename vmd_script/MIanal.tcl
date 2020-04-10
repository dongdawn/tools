proc MIanal {outfile molid nHIS HISframe Bprcnt blength selPtxt selNtxt} {
### option ###
### 其中，hbABnum为A-H-B氢键总数，hbListBA为B-H-A氢键总数
###				HBList($unt)为$unt原子涉及氢键的次数
###				HBFrmList($unt$chF$iHIS)为第iHIS时间段内$unt原子涉及氢键的次数
###				HBpairList($hbH$chA$hbA)为某氢键出现的总次数
###				HBpairFrmList($hbH$chA$hbA$chF$iHIS)为第iHIS时间段内某氢键出现的次数
  set outFL [open $outfile w]
  set hz HIS
  set outHISFL [open $outfile$hz w]
  set pr pair
  set outHISprFL [open $outfile$hz$pr w]
  set outHISprPFL [open $outfile$hz${pr}P w]
  set cutoff 30.0
#  set molid [molinfo top]
  set selP [atomselect $molid $selPtxt]
  set selN [atomselect $molid $selNtxt]
  set selALL [atomselect $molid all]
  set pairP ""
  set pairN "" 
  foreach uip [$selP get index] uchp [$selP get chain] uidp [$selP get resid] {
    foreach uin [$selN get index] uchn [$selN get chain] uidn [$selN get resid] {
      if {$uchp != $uchn || $uidp != $uidn} then {
      if {[measure bond "$uip $uin" molid $molid frame 0] < $cutoff} then {
        lappend pairP $uip
        lappend pairN $uin
      }
      }
    }
  }
  set nframe [molinfo $molid get numframes]
#  set HISframe [expr "$nframe/$nHIS"]
  set Blow [expr "int($nframe*$Bprcnt)"]
  set chA A
  set chF F
  if {[info exists BFrmList]} {unset BFrmList}
  if {[info exists BpairList]} {unset BpairList}
  if {[info exists BpairFrmList]} {unset BpairFrmList}
## loop over frames ##
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    set frmsum 0
    foreach up $pairP un $pairN {
      if {[measure bond "$up $un" molid $molid frame $iframe] <= $blength} then {
        incr frmsum
        set iHIS [expr "int($iframe/$HISframe)"]
        if {![info exists BFrmList($up$chF$iHIS)]} {set BFrmList($up$chF$iHIS) 0}
        incr BFrmList($up$chF$iHIS)
        if {![info exists BFrmList($un$chF$iHIS)]} {set BFrmList($un$chF$iHIS) 0}
        incr BFrmList($un$chF$iHIS)
        if {![info exists BpairList($up$chA$un)]} {set BpairList($up$chA$un) 0}
        incr BpairList($up$chA$un)
        if {![info exists BpairFrmList($up$chA$un$chF$iHIS)]} {set BpairFrmList($up$chA$un$chF$iHIS) 0}
        incr BpairFrmList($up$chA$un$chF$iHIS)
      }
    }
    puts $outFL "$iframe $frmsum"
  }
## puts BFrmList ##
  foreach atmU [lsort -unique -integer "$pairP $pairN"] {
    set outtxt ""
    lappend outtxt $atmU
    lappend outtxt [lindex [$selALL get chain] $atmU]
    lappend outtxt [lindex [$selALL get resname] $atmU]
    lappend outtxt [lindex [$selALL get resid] $atmU]
    lappend outtxt [lindex [$selALL get name] $atmU]
    lappend outtxt :
    for {set iHIS 0} {$iHIS < $nHIS} {incr iHIS} {
      if {![info exists BFrmList($atmU$chF$iHIS)]} {set BFrmList($atmU$chF$iHIS) 0}
      lappend outtxt $BFrmList($atmU$chF$iHIS)
    }
    puts $outHISFL $outtxt
  }
# puts BpairList
  foreach up $pairP un $pairN {
    set outtxt ""
    if {[info exists BpairList($up$chA$un)] && $BpairList($up$chA$un) > $Blow} {
      lappend outtxt $up
      lappend outtxt [lindex [$selALL get chain] $up]
      lappend outtxt [lindex [$selALL get resname] $up]
      lappend outtxt [lindex [$selALL get resid] $up]
      lappend outtxt [lindex [$selALL get name] $up]
      lappend outtxt -
      lappend outtxt $un
      lappend outtxt [lindex [$selALL get chain] $un]
      lappend outtxt [lindex [$selALL get resname] $un]
      lappend outtxt [lindex [$selALL get resid] $un]
      lappend outtxt [lindex [$selALL get name] $un]
      set outtxtp $outtxt 
      puts $outHISprPFL "$outtxtp [expr $BpairList($up$chA$un) / double($nframe)]"
      lappend outtxt $BpairList($up$chA$un)
      lappend outtxt :
      for {set iHIS 0} {$iHIS < $nHIS} {incr iHIS} {
        if {![info exists BpairFrmList($up$chA$un$chF$iHIS)]} {set BpairFrmList($up$chA$un$chF$iHIS) 0}
        lappend outtxt $BpairFrmList($up$chA$un$chF$iHIS)
      }
      puts $outHISprFL $outtxt
    }
  }
  close $outFL
  close $outHISFL
  close $outHISprFL
  close $outHISprPFL
}
#MIanal stBG top 50.0 200.0 0.1 4.0 "(resname ARG and name NH1 NH2) or (resname LYS and name NZ) or (resname HSP and name ND1 NE2) or (name N and (same residue as name HT3 H3))" "(resname GLU and name OE2 OE1) or (resname ASP and name OD2 OD1) or name OT1 OT2 O1 O2"
#MIanal stBGdpcN top 50.0 200.0 0.1 4.0 "resname DPC and name N" "(resname GLU and name OE2 OE1) or (resname ASP and name OD2 OD1) or name OT1 OT2 O1 O2"
#MIanal stBGdpcP top 50.0 200.0 0.1 4.0 "(resname ARG and name NH1 NH2) or (resname LYS and name NZ) or (resname HSP and name ND1 NE2) or (name N and (same residue as name HT3 H3))" "resname DPC and name O31 O32 O33 O34"