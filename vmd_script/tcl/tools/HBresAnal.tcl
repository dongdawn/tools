# 20111223 V2.2: faster by defining more and more detailed selections #
# 20111223 V2.1: faster by defining more detailed selections #
#                sel1 and sel2 must IDENTICAL or ABSULUTELY DIFFERENT, or there will be error caused by VMD
# 20111221 V2.0: seltxt1 & seltxt2 can be the same #
proc HBresAnal {outfile molid HBprcnt HBlength HBangle seltxt1 seltxt2} {
### option ###
### 所有帧的分段个数nHIS，出现频率底限HBprcnt ###
### 其中，hbABnum为A-H-B氢键总数，hbListBA为B-H-A氢键总数
###				HBpairList($hbH$cA$hbA)为某氢键出现的总次数
###				HBpairFrmList($hbH$cA$hbA$cF$iHIS)为第iHIS时间段内某氢键出现的次数
## parameters ##
  set wfl [open ${outfile}RESP w]
  set nframe [molinfo $molid get numframes]
  set HBlow [expr "int($nframe*$HBprcnt)"]
  set selALL [atomselect $molid all]
  set cA A
  set cF F
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
    foreach urs $resSEList uds $dmnSEList {
      $urs frame $iframe
      $uds frame $iframe
      set hbListAB [measure hbonds $HBlength $HBangle $urs $uds]
      set hbListBA [measure hbonds $HBlength $HBangle $uds $urs]
      foreach hbA [lindex $hbListAB 1] hbH [lindex $hbListAB 2] {
        if {![info exists HBpairList($hbH$cA$hbA)]} {set HBpairList($hbH$cA$hbA) 0}
        incr HBpairList($hbH$cA$hbA)
        if {![info exists HBpairFrmList($hbH$cA$hbA$cF$iframe)]} {set HBpairFrmList($hbH$cA$hbA$cF$iframe) 0}
        incr HBpairFrmList($hbH$cA$hbA$cF$iframe)
      }
      foreach hbA [lindex $hbListBA 1] hbH [lindex $hbListBA 2] {
        if {![info exists HBpairList($hbA$cA$hbH)]} {set HBpairList($hbA$cA$hbH) 0}
        incr HBpairList($hbA$cA$hbH)
        if {![info exists HBpairFrmList($hbA$cA$hbH$cF$iframe)]} {set HBpairFrmList($hbA$cA$hbH$cF$iframe) 0}
        incr HBpairFrmList($hbA$cA$hbH$cF$iframe)
      }
    }
  }
  if {[info exists RESpair1]} {unset RESpair1}
  if {[info exists RESpair2]} {unset RESpair2}
  if {[info exists RESpairFn]} {unset RESpairFn}
  set RESpairs ""
  foreach {nm vl} [array get HBpairList] {
    set whA [string first $cA $nm]
    set nm1 [string range $nm 0 [expr $whA-1]]
    set nm2 [string range $nm [expr $whA+1] [expr [string length $nm]-1]]
    if {$vl > $HBlow} {
      set tp "[lindex [$selALL get chain] $nm1][lindex [$selALL get resid] $nm1][lindex [$selALL get chain] $nm2][lindex [$selALL get resid] $nm2]"
      set RESpairs "$RESpairs $tp"
      set RESpair1($tp) $nm1
      set RESpair2($tp) $nm2
      for {set iframe 0} {$iframe < $nframe} {incr iframe} {
        if {[info exists HBpairFrmList($nm1$cA$nm2$cF$iframe)] && $HBpairFrmList($nm1$cA$nm2$cF$iframe)==1} then {
          set RESpairFn($tp,$iframe) 1
        }
      }
    }
  }
  set RESpairs [lsort -unique $RESpairs]
  foreach ur $RESpairs {
    for {set iframe 0} {$iframe < $nframe} {incr iframe} {
      if {![info exists RESpairFn($ur,$iframe)]} {set RESpairFn($ur,$iframe) 0}
    }
  }
  foreach ur $RESpairs {
    set RESpairTF($ur) 0
    for {set iframe 0} {$iframe < $nframe} {incr iframe} {
      if {$RESpairFn($ur,$iframe)==1} {incr RESpairTF($ur)}
    }
  }
# puts RESpairs
  foreach ur $RESpairs {
    set outxt ""
    lappend outxt [lindex [$selALL get chain] $RESpair1($ur)]
    lappend outxt [lindex [$selALL get resname] $RESpair1($ur)]
    lappend outxt [lindex [$selALL get resid] $RESpair1($ur)]
    lappend outxt -
    lappend outxt [lindex [$selALL get chain] $RESpair2($ur)]
    lappend outxt [lindex [$selALL get resname] $RESpair2($ur)]
    lappend outxt [lindex [$selALL get resid] $RESpair2($ur)]
    puts $wfl "$outxt [expr double($RESpairTF($ur)) / double($nframe)]"
  }
  close $wfl
}
# outfile molid HBprcnt HBlength HBangle seltxt1 seltxt2
# HBresAnal chAB-FG $molid 0.0 3.5 30.0 "chain A B" "chain F G"
#HBresAnal C2_275-308 0 0.0 3.5 30.0 "chain C D and resid 275" "chain C D and resid 308"
#HBresAnal C2_276-308 0 0.0 3.5 30.0 "chain C D and resid 276" "chain C D and resid 308"
#HBresAnal C2_279-308 0 0.0 3.5 30.0 "chain C D and resid 279" "chain C D and resid 308"
#HBresAnal C2_278-306 0 0.0 3.5 30.0 "chain C D and resid 306" "chain C D and resid 278"
#HBresAnal C2_252-303 0 0.0 3.5 30.0 "chain C D and resid 252" "chain C D and resid 303"
#HBresAnal C2_252-305 0 0.0 3.5 30.0 "chain C D and resid 252" "chain C D and resid 305"
#HBresAnal C2_252-306 0 0.0 3.5 30.0 "chain C D and resid 252" "chain C D and resid 306"
#HBresAnal C2_253-279 0 0.0 3.5 30.0 "chain C D and resid 253" "chain C D and resid 279"
#HBresAnal C2_205-302 0 0.0 3.5 30.0 "chain C D and resid 205" "chain C D and resid 302"
#HBresAnal resid_500 0 0.0 3.5 30.0 "resid 500" "protein"


