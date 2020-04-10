# cd {F:\TolC\MSM\pH7.5\MSM\round2}
set flname rmsd-mutual
set seltxt "name CA"
set pdbfile {E:\GGBP\4.after\2fw0_start.pdb}
set trjpath {E:\GGBP\4.after}
set maxrmsd 11.
set minrmsd 0.
set drmsd 0.1
set step 20

set utrjid [mol new $pdbfile type pdb]
set vtrjid [mol new $pdbfile type pdb]
set usel [atomselect $utrjid $seltxt]
set vsel [atomselect $vtrjid $seltxt]
animate delete beg 0 end 0 $utrjid
animate delete beg 0 end 0 $vtrjid

if {[info exists count] == 1} {unset count}
if {[info exists wline] == 1} {unset wline}
for {set u $minrmsd} {$u<=$maxrmsd} {set u [expr $u + $drmsd]} {
  set tp [expr round($u * 10)/ 10.0]
  set wline($tp) $tp
  set count($tp) 0
}

cd $trjpath
#[glob *]
set utrajlist {2fvy01 2fvy02 2fw001 2fw002 2fw003 2gbp01 2gbp02 2gbp03 2hph01 2hph02 2hph03 2qw101 2qw102 2qw103 }
set vtrajlist {2fvy01 2fvy02 2fw001 2fw002 2fw003 2gbp01 2gbp02 2gbp03 2hph01 2hph02 2hph03 2qw101 2qw102 2qw103 }
#set utrajlist {WT_f}
#set vtrajlist {D153A_14_1a D153A_14_1b D153A_14_1c D153A_14_2a D153A_14_2b D153A_14_2c D153A_14_3a D153A_14_3b D153A_14_3c D153A_9_1a D153A_9_1b D153A_9_2a D153A_9_2b D153A_9_3a D153A_9_3b TolC0 TolC1 TolC10 TolC11 TolC12 TolC13 TolC14 TolC15 TolC16 TolC17 TolC18 TolC19 TolC2 TolC20 TolC21 TolC22 TolC23 TolC24 TolC25 TolC26 TolC27 TolC28 TolC29 TolC3 TolC30 TolC31 TolC32 TolC33 TolC34 TolC35 TolC36 TolC37 TolC38 TolC39 TolC4 TolC40 TolC41 TolC42 TolC43 TolC44 TolC45 TolC46 TolC47 TolC48 TolC49 TolC5 TolC6 TolC7 TolC8 TolC9 WT_2_1a WT_2_1b WT_2_1c WT_2_2a WT_2_2b WT_2_2c WT_2_3a WT_2_3b WT_2_3c WT_3_1a WT_3_1b WT_3_1c WT_3_1d WT_3_2a WT_3_2b WT_3_2c WT_3_2d WT_3_3a WT_3_3b WT_3_3c WT_3_3d Y362FR367S_6_1a Y362FR367S_6_1b Y362FR367S_6_1c Y362FR367S_6_2a Y362FR367S_6_2b Y362FR367S_6_2c Y362FR367S_6_3a Y362FR367S_6_3b Y362FR367S_6_3c}
#set utrajlist {WT_c}
#set vtrajlist {WT_c}
foreach ut $utrajlist {
  set fl [open ${flname}_$ut w]
  cd $ut
  mol top $utrjid
  mol addfile md_noPBC.xtc type xtc first 5000 last -1 step $step filebonds 1 autobonds 1 waitfor all
  set utrjsel [atomselect $utrjid $seltxt]
  set unframe [molinfo $utrjid get numframes]
  cd ..
  set uheadline bin
  set vheadline bin
  foreach vt $vtrajlist {
    lappend uheadline $ut
    lappend vheadline $vt
    for {set u $minrmsd} {$u<=$maxrmsd} {set u [expr $u + $drmsd]} {
      set tp [expr round($u * 10)/ 10.0]
      set count($tp) 0
    }
    cd $vt
    mol top $vtrjid
    mol addfile md_noPBC.xtc type xtc first 0 last -1 step $step filebonds 1 autobonds 1 waitfor all
    set vtrjsel [atomselect $vtrjid $seltxt]
    set vnframe [molinfo $vtrjid get numframes]
    for {set iframe 0} {$iframe < $unframe} {incr iframe} {
      for {set jframe 0} {$jframe < $vnframe} {incr jframe} {
        $utrjsel frame $iframe
        $vtrjsel frame $jframe
        set trans_mat [measure fit $utrjsel $vtrjsel weight mass]
        $utrjsel move $trans_mat
        set rmsd [measure rmsd $utrjsel $vtrjsel weight mass]
        set tp [expr round($rmsd * 10) / 10.0]
        incr count($tp)
      }
    }
    for {set u $minrmsd} {$u<=$maxrmsd} {set u [expr $u + $drmsd]} {
      set tp [expr round($u * 10)/ 10.0]
      lappend wline($tp) $count($tp)
    }
    animate delete beg 0 end [expr $vnframe -1] $vtrjid
    cd ..
  }
  animate delete beg 0 end [expr $unframe -1] $utrjid
  puts $fl $uheadline
  puts $fl $vheadline
  for {set u $minrmsd} {$u<=$maxrmsd} {set u [expr $u + $drmsd]} {
    set tp [expr round($u * 10)/ 10.0]
    puts $fl $wline($tp)
    set wline($tp) $tp
    set count($tp) 0
  }
  close $fl
}


