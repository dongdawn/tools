proc MIresAnal {outfile molid Bprcnt blength selPtxt selNtxt} {
### option ###
### 其中，BList($unt)为$unt原子涉及氢键的次数
###				BpairList($hbH$chA$hbA)为某氢键出现的总次数
###				BpairFrmList($hbH$chA$hbA$chF$iHIS)为第iHIS时间段内某氢键出现的次数
  set wfl [open ${outfile}RESP w]
  set cutoff 200.0
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
  set Blow [expr "int($nframe*$Bprcnt)"]
  set chA A
  set chF F
  if {[info exists BpairList]} {unset BpairList}
  if {[info exists BpairFrmList]} {unset BpairFrmList}
## loop over frames ##
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    foreach up $pairP un $pairN {
      if {[measure bond "$up $un" molid $molid frame $iframe] <= $blength} then {
        if {![info exists BpairList($up$chA$un)]} {set BpairList($up$chA$un) 0}
        incr BpairList($up$chA$un)
        if {![info exists BpairFrmList($up$chA$un$chF$iframe)]} {set BpairFrmList($up$chA$un$chF$iframe) 0}
        incr BpairFrmList($up$chA$un$chF$iframe)
      }
    }
  }
  if {[info exists RESpair1]} {unset RESpair1}
  if {[info exists RESpair2]} {unset RESpair2}
  if {[info exists RESpairFn]} {unset RESpairFn}
  set RESpairs ""
  foreach up $pairP un $pairN {
    if {[info exists BpairList($up$chA$un)] && $BpairList($up$chA$un) > $Blow} {
      if {[lindex [$selALL get chain] $up] < [lindex [$selALL get chain] $un]} then {
        set tp "[lindex [$selALL get chain] $up][lindex [$selALL get resid] $up][lindex [$selALL get chain] $un][lindex [$selALL get resid] $un]"
      } else {
        set tp "[lindex [$selALL get chain] $un][lindex [$selALL get resid] $un][lindex [$selALL get chain] $up][lindex [$selALL get resid] $up]"
      }
      set RESpairs "$RESpairs $tp"
      set RESpair1($tp) $up
      set RESpair2($tp) $un
      for {set iframe 0} {$iframe < $nframe} {incr iframe} {
        if {[info exists BpairFrmList($up$chA$un$chF$iframe)] && $BpairFrmList($up$chA$un$chF$iframe)==1} then {
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
    #puts $wfl "$outxt [expr double($RESpairTF($ur)) / double($nframe)]"
    puts $wfl "$outxt [expr double($RESpairTF($ur)) / double($nframe)]     expr abs([lindex [$selALL get resid] $RESpair1($ur)]-[lindex [$selALL get resid] $RESpair2($ur)])     $RESpairTF($ur)"
  }
  close $wfl
}
#MIresAnal stBG_GCFC $molid 0.0 4.0 "(resname ARG and name NH1 NH2) or (resname LYS and name NZ) or (resname HSP and name ND1 NE2) or (name N and (same residue as name HT3))" "(chain F and resid 514 and name OD1 OD2 OT1 OT2) or (chain G and resid 296 and name OT1 OT2)"
#MIresAnal met1hph 0 0.0 4.0 "resid 345 and name CB CG SD CE" "not backbone and noh and resid 273 278 282 296 301 312"
#MIresAnal met2hph 0 0.0 4.0 "resid 346 and name CB CG SD CE" "not backbone and noh and resid 273 278 282 296 301 312"
#MIresAnal met301hph 0 0.0 4.0 "resid 301 and name CB CG SD CE" "not backbone and noh and not name \"N.*\" \"O.*\" and not resid 301 and not water"
#MIresAnal met312hph 0 0.0 4.0 "resid 312 and name CB CG SD CE" "not backbone and noh and not name \"N.*\" \"O.*\" and not resid 312 and not water"
#MIresAnal K308_6 0 0.0 6.0 "resid 308 and name NZ" "noh and not resid 308 and not water"
#MIresAnal K308_4 0 0.0 4.0 "resid 308 and name NZ" "noh and not resid 308 and not water"
#MIresAnal K308_6bb 0 0.0 6.0 "resid 308 and name NZ" "noh and not resid 308 and not water and backbone"
#MIresAnal K308_4bb 0 0.0 4.0 "resid 308 and name NZ" "noh and not resid 308 and not water and backbone"
#MIresAnal K308_45bb 0 0.0 4.5 "resid 308 and name NZ" "noh and not resid 308 and not water and backbone"
#MIresAnal K308_5bb  0 0.0 5.0 "resid 308 and name NZ" "noh and not resid 308 and not water and backbone"
   