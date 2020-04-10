# do not move REF before using this
# create new PDB file and remove it after used, so none loaded files are moved
proc rmsdDMN {outfile molidREF molidTRAJ fstep sel} {
  set refSNAPname rmsdDMNtmp
  set refSET [atomselect $molidREF "$sel" frame 0]
  if {[$refSET num]!=0} then {
    set fl [open $outfile w]
    $refSET writepdb $refSNAPname
    set refID [mol load pdb $refSNAPname]
    set ref [atomselect $refID "$sel" frame 0]
    set traj [atomselect $molidTRAJ "$sel"]
    set nframe [molinfo $molidTRAJ get numframes]
    for {set iframe 0} {$iframe < $nframe} {incr iframe $fstep} {
      $traj frame $iframe
      set trans_mat [measure fit $ref $traj weight mass]
      $ref move $trans_mat
      puts $fl "$iframe [measure rmsd $traj $ref weight mass]"
    }
    mol delete $refID
    mol top $molidTRAJ
    close $fl
  }
#  del rmsdDMNtmp
}
#rmsdDMN rmsd2ptNBDz 0 top 1 "name CA and chain A B and resid 1 to 235"
#rmsdDMN rmsd2ptFGcorez 0 top 1 "((chain G and resid 80 to 283) or (chain F and resid 276 to 505)) and name CA"

# align traj to REF and calc part of atoms in REF
proc rmsdPart {outfile molidREF molidTRAJ selOL selCalc} {
  set refSNAPname rmsdPtmp
  set fl [open $outfile w]
  [atomselect $molidREF all frame 0] writepdb $refSNAPname
  set refID [mol load pdb ./$refSNAPname]
  mol top 0
  set ref [atomselect $refID "$selOL" frame 0]
  set traj [atomselect $molidTRAJ "$selOL"]
  set refCalc [atomselect $refID "$selCalc" frame 0]
  set trajCalc [atomselect $molidTRAJ "$selCalc"]
  set refAll [atomselect $refID all frame 0]
  set nframe [molinfo $molidTRAJ get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $traj frame $iframe
    $trajCalc frame $iframe
    set trans_mat [measure fit $ref $traj weight mass]
    $refAll move $trans_mat
    puts $fl "$iframe [measure rmsd $trajCalc $refCalc weight mass]"
  }
  mol delete $refID
  close $fl
}
#rmsdPart rmsd2initCA_tm345A 1 0 "name CA" "name CA and chain A and resid 92 to 168"
#rmsdPart rmsd2initCA_tm345B 1 0 "name CA" "name CA and chain B and resid 92 to 168"

# move TRAJ
proc alignTRAJ {refID trajID sel} {
  set ref [atomselect $refID "$sel" frame 0]
  set traj [atomselect $trajID "$sel"]
  set trajALL [atomselect $trajID all]
  set nframe [molinfo $trajID get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $traj frame $iframe
    $trajALL frame $iframe
    set trans_mat [measure fit $traj $ref weight mass]
    $trajALL move $trans_mat
  }
}
#alignTRAJ 1 0 "name CA"

proc alignPmovQ {refID trajID fitSEL movSEL} {
  set ref [atomselect $refID $fitSEL frame 0]
  set traj [atomselect $trajID $fitSEL]
  set movtraj [atomselect $trajID $movSEL]
  set nframe [molinfo $trajID get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $traj frame $iframe
    $movtraj frame $iframe
    set trans_mat [measure fit $traj $ref weight mass]
    $movtraj move $trans_mat
  }
}
#alignTRAJ 1 0 "name CA" "all"

# align trajALL by fitting selTRAJ to selREF
proc alignSEL {refID selREF trajID selTRAJ} {
  set ref [atomselect $refID "$selREF" frame 0]
  set traj [atomselect $trajID "$selTRAJ"]
  set trajALL [atomselect $trajID all]
  set nframe [molinfo $trajID get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $traj frame $iframe
    $trajALL frame $iframe
    set trans_mat [measure fit $traj $ref weight mass]
    $trajALL move $trans_mat
  }
}
# PDZ4-PDZ5
#alignSEL 1 "name CA and resid 485 to 494 495 to 507 513 to 525 528 529 530 to 548 551 to 566 568 to 579" \
#         0 "name CA and resid 580 to 589 591 to 603 604 to 616 617 618 621 to 639 640 to 655 656 to 667"



