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

#       输出文件名  分子ID1           选区一              分子ID2            选取二
MinDist mdist_64-184 0 "chain A and resid 64" 0 "chain A and resid 184"
MinDist mdist_64-185 0 "chain A and resid 64" 0 "chain A and resid 185"
MinDist mdist_39-157 0 "chain A and resid 39" 0 "chain A and resid 157"
MinDist mdist_43-66  0 "chain A and resid 43" 0 "chain A and resid 66"
MinDist mdist_40-184 0 "chain A and resid 40" 0 "chain A and resid 184"
MinDist mdist_39-167 0 "chain A and resid 39" 0 "chain A and resid 167"

MCenterDist distCA_39-157 0 "chain A and resid 39 and name CA" "chain A and resid 157 and name CA"
MCenterDist distCA_39-167 0 "chain A and resid 39 and name CA" "chain A and resid 167 and name CA"
MCenterDist distCA_40-184 0 "chain A and resid 40 and name CA" "chain A and resid 184 and name CA"
MCenterDist distCA_40-185 0 "chain A and resid 40 and name CA" "chain A and resid 185 and name CA"
MCenterDist distCA_40-64  0 "chain A and resid 40 and name CA" "chain A and resid 64 and name CA"
MCenterDist distCA_44-66  0 "chain A and resid 44 and name CA" "chain A and resid 66 and name CA"

rmsdDMN rmsd2i  1  0  1 "chain A and name CA"
rmsdDMN rmsd2isnr 1 0 1 "chain A and name CA and resid 138 to 195"
rmsdDMN rmsd2ilgn 1 0 1 "chain A and name CA and resid 1 to 137"


set selcalcList {"chain A and name CA" "chain A and name CA and resid 1 to 137" "chain A and name CA and resid 138 to 195"}
set selovaList  {"chain A and name CA" "chain A and name CA and resid 1 to 137" "chain A and name CA and resid 138 to 195"}
set flnmList    { rmsfALL               rmsfLGN                                  rmsfSNR                                  }
#         selcalcList  selovaList  flnmList refID molID initF lastF
calcrmsf $selcalcList $selovaList $flnmList     1     0     0 44999


proc ATMdihed {outfile molid fstep sel1 sel2 sel3 sel4} {
  set fl [open $outfile w]
  set motif1 [atomselect $molid "$sel1"]
  set motif2 [atomselect $molid "$sel2"]
  set motif3 [atomselect $molid "$sel3"]
  set motif4 [atomselect $molid "$sel4"]
  if {[$motif1 num]==1 && [$motif2 num]==1 && [$motif3 num]==1 && [$motif4 num]==1} then {
    set atm1 [$motif1 get index]
    set atm2 [$motif2 get index]
    set atm3 [$motif3 get index]
    set atm4 [$motif4 get index]
    set nframe [molinfo $molid get numframes]
    for {set iframe 0} {$iframe < $nframe} {incr iframe $fstep} {
      set tp [measure dihed "$atm1 $atm2 $atm3 $atm4" molid $molid frame $iframe]
      puts $fl "$iframe $tp"
    }
  } else {
    puts "Only one atom per selection is allowed!"
  }
  close $fl
}

ATMdihe dihed4w 0 1 "chain A and resid 4 and name N" "chain A and resid 4 and name CA" "chain A and resid 4 and name CB" "chain A and resid 4 and name CG"
ATMdihe dihed4x 0 1 "chain A and resid 4 and name CA" "chain A and resid 4 and name CB" "chain A and resid 4 and name CG" "chain A and resid 4 and name CD1"

ATMdihe dihed28w 0 1 "chain A and resid 28 and name N" "chain A and resid 28 and name CA" "chain A and resid 28 and name CB" "chain A and resid 28 and name CG"
ATMdihe dihed28x 0 1 "chain A and resid 28 and name CA" "chain A and resid 28 and name CB" "chain A and resid 28 and name CG" "chain A and resid 28 and name CD1"

ATMdihe dihed31w 0 1 "chain A and resid 31 and name N" "chain A and resid 31 and name CA" "chain A and resid 31 and name CB" "chain A and resid 31 and name CG"
ATMdihe dihed31x 0 1 "chain A and resid 31 and name CA" "chain A and resid 31 and name CB" "chain A and resid 31 and name CG" "chain A and resid 31 and name CD1"

ATMdihe dihed39w 0 1 "chain A and resid 39 and name N" "chain A and resid 39 and name CA" "chain A and resid 39 and name CB" "chain A and resid 39 and name CG"
ATMdihe dihed39x 0 1 "chain A and resid 39 and name CA" "chain A and resid 39 and name CB" "chain A and resid 39 and name CG" "chain A and resid 39 and name CD1"

ATMdihe dihed42w 0 1 "chain A and resid 42 and name N" "chain A and resid 42 and name CA" "chain A and resid 42 and name CB" "chain A and resid 42 and name CG"
ATMdihe dihed42x 0 1 "chain A and resid 42 and name CA" "chain A and resid 42 and name CB" "chain A and resid 42 and name CG" "chain A and resid 42 and name CD1"

ATMdihe dihed64w 0 1 "chain A and resid 64 and name N" "chain A and resid 64 and name CA" "chain A and resid 64 and name CB" "chain A and resid 64 and name CG"
ATMdihe dihed64x 0 1 "chain A and resid 64 and name CA" "chain A and resid 64 and name CB" "chain A and resid 64 and name CG" "chain A and resid 64 and name CD1"

ATMdihe dihed184w 0 1 "chain A and resid 184 and name N" "chain A and resid 184 and name CA" "chain A and resid 184 and name CB" "chain A and resid 184 and name CG"
ATMdihe dihed184x 0 1 "chain A and resid 184 and name CA" "chain A and resid 184 and name CB" "chain A and resid 184 and name CG" "chain A and resid 184 and name CD1" 

ATMdihe dihed185w 0 1 "chain A and resid 185 and name N" "chain A and resid 185 and name CA" "chain A and resid 185 and name CB" "chain A and resid 185 and name CG"
ATMdihe dihed185x 0 1 "chain A and resid 185 and name CA" "chain A and resid 185 and name CB" "chain A and resid 185 and name CG" "chain A and resid 185 and name CD1" 







proc calcBSA {filename molID solrad initF lastF sel1txt sel2txt} {
# parameters #
  set fl [open $filename w]
  set sel12 [atomselect $molID "($sel1txt) or ($sel2txt)"]
  set sel1 [atomselect $molID $sel1txt]
  set sel2 [atomselect $molID $sel2txt]
  set nframe [molinfo $molID get numframes]
  if {$lastF == -1} then { set lastF [expr $nframe -1] }
# loop in frames #
  for {set iframe $initF} {$iframe <= $lastF} {incr iframe 10} {
    $sel12 frame $iframe
    $sel1 frame $iframe
    $sel2 frame $iframe
    set sa12 [measure sasa $solrad $sel12]
    set sa1 [measure sasa $solrad $sel1]
    set sa2 [measure sasa $solrad $sel2]
    puts $fl "$iframe [expr $sa1 + $sa2 - $sa12] $sa1 $sa2 $sa12"
  }
  close $fl 
}
#
calcBSA BSAlgn-snaredpc 0 1.4 0 8999 "resid 1 to 137" "resid 138 to 195 or resname DPC"
calcBSA BSAlgn-snarefar 0 1.4 0 8999 "resid 1 to 137" "resid 138 to 195 or resname FAR"
calcBSA BSAlgn-dpc 0 1.4 0 8999 "resid 1 to 137" "resname DPC"
calcBSA BSAlgn-far 0 1.4 0 8999 "resid 1 to 137" "resname FAR"
calcBSA BSAlgn-aFG   0 1.4 0 8999 "resid 1 to 137" "resid 167 to 192"
calcBSA BSAlgn-aDEBF   0 1.4 0 8999 "resid 1 to 137" "resid 141 to 160"
calcBSA BSAlgn-snare 0 1.4 0 8999 "resid 1 to 137" "resid 138 to 195"
calcBSA BSAlgn-snareprot 0 1.4 0 8999 "resid 1 to 137" "resid 138 to 194"
calcBSA BSAlgn-aE    0 1.4 0 8999 "resid 1 to 137" "resid 157 to 167"
calcBSA BSAlgn-BF    0 1.4 0 8999 "resid 1 to 137" "resid 152 to 156"
calcBSA BSAlgn-aD    0 1.4 0 8999 "resid 1 to 137" "resid 141 to 151"






# based on STRIDE used by VMD 
# OUTPUT: 1.ss_percent: ss percentage of the whole protein in each bin
# OUTPUT: 2.ss_percentRES: ss percentage of each residue
# helix=H+I+G
# sheet=E+B
# T
# C

proc ssappear {list} {
  set ssH 0
  set ssI 0
  set ssG 0
  set ssE 0
  set ssB 0
  set ssT 0
  set ssC 0
  foreach u $list {
    switch $u {
      H {incr ssH}
      I {incr ssI}
      G {incr ssG}
      E {incr ssE}
      B {incr ssB}
      T {incr ssT}
      C {incr ssC}
    }
  }
  return "$ssH $ssI $ssG $ssE $ssB $ssT $ssC"
}

### options ###
set trajID 0
set frminbin 100
set outfile ss_percent
set outfileres ss_percentRES
set seltxt "name CA"

### parameters ###
set ssName {     H         I         G            E            B       T    C  }
set ssList {alpha_helix pi_helix helix_3_10 extended_beta bridge_beta turn coil}
set traj [atomselect $trajID $seltxt]
set nres [$traj num]
set fl [open $outfile w]
set flrs [open $outfileres w]

### initialize ###
if {[info exists bin]} then {unset bin}
set nframe [molinfo $trajID get numframes]
set ibin 0
set ifrminbin 0
foreach sn $ssName {set bin(${sn}_${ibin}) 0}
set ires 0
foreach uch [$traj get chain] urn [$traj get resname] uri [$traj get resid] {
  set rinfo($ires) "$uch $urn $uri"
  set ss_res($ires) ""
  incr ires
}

### loop in frame ###
for {set iframe 0} {$iframe < $nframe} {incr iframe} {
  incr ifrminbin
  if {$ifrminbin > $frminbin} then {
    set fib($ibin) $ifrminbin
    incr ibin
    foreach sn $ssName {set bin(${sn}_${ibin}) 0}
    set ifrminbin 1
  }
  molinfo $trajID set frame $iframe
  mol ssrecalc $trajID
  set sstrj [$traj get structure]
  foreach sn $ssName snum [ssappear $sstrj] {
    set bin(${sn}_${ibin}) [expr $bin(${sn}_${ibin}) + $snum]
  }
  set ires 0
  foreach us $sstrj {
    lappend ss_res($ires) $us
    incr ires
  }
}
set fib($ibin) $ifrminbin

### output ###
set nbin $ibin
puts $fl "bin # frame # $ssList"
for {set ibin 0} {$ibin <= $nbin} {incr ibin} {
  set outtxt ""
  foreach sn $ssName {
    lappend outtxt [expr 1.0 * $bin(${sn}_${ibin}) / ( $bin(H_${ibin}) + $bin(I_${ibin}) + $bin(G_${ibin}) + $bin(E_${ibin}) + $bin(B_${ibin}) + $bin(T_${ibin}) + $bin(C_${ibin}) ) ]
  }
  puts $fl "$ibin [expr $ibin * $frminbin] - [expr $ibin * $frminbin + $fib($ibin) -1] $outtxt"
}
close $fl

puts $flrs "chn rnm rid $ssList"
for {set ires 0} {$ires < $nres} {incr ires} {
  set outxt $rinfo($ires)
  foreach us [ssappear $ss_res($ires)] {
    set outxt "$outxt [expr $us / double($nframe)]"
  }
  puts $flrs $outxt
}
close $flrs

