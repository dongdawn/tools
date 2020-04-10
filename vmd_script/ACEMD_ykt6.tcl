set tclList {"E:/weng/tcl/multiDCD/tools/Ax-Ax.tcl" "E:/weng/tcl/multiDCD/tools/HBpairAnal.tcl" "E:/weng/tcl/multiDCD/tools/MCdist.tcl" "E:/weng/tcl/multiDCD/tools/rmsd.tcl" "E:/weng/tcl/multiDCD/tools/MIanal.tcl" "E:/weng/tcl/multiDCD/tools/MIresAnal.tcl" "E:/weng/tcl/multiDCD/tools/HBresAnal.tcl" "E:/weng/tcl/multiDCD/tools/rmsf.tcl"}

foreach ut $tclList {
  source $ut
}

set molid 0

rmsdDMN rmsd2i    0 $molid 1 "chain A and name CA"
rmsdDMN rmsd2ilgn 0 $molid 1 "chain A and name CA and resid 1 to 137"
rmsdDMN rmsd2isnr 0 $molid 1 "chain A and name CA and resid 138 to 195"

MCenterDist dist_K60hn-H154o $molid "chain A and resid 60 and name HN" "chain A and resid 154 and name O"
MCenterDist dist_S58o-H154hn $molid "chain A and resid 58 and name O"  "chain A and resid 154 and name HN"
MCenterDist dist_S58hn-H152o $molid "chain A and resid 58 and name HN" "chain A and resid 152 and name O"
MCenterDist dist_R56o-I152hn $molid "chain A and resid 56 and name O"  "chain A and resid 152 and name HN"

set selcalcList {"chain A and name CA" "chain A and name CA and resid 1 to 137" "chain A and name CA and resid 138 to 195"}
set selovaList  {"chain A and name CA" "chain A and name CA and resid 1 to 137" "chain A and name CA and resid 138 to 195"}
set flnmList    { rmsfALL               rmsfLGN                                  rmsfSNR                                  }
#         selcalcList  selovaList  flnmList refID molID initF lastF
calcrmsf $selcalcList $selovaList $flnmList     1     0     0  8999

# 9000֡ prot10.vmd 1.8 us
HBpairAnal lgn-snr 100.0 90.0 0.05 $molid 3.5 30.0 "chain A and resid 1 to 137" "chain A and resid 138 to 195"
HBpairAnal lgn 100.0 90.0 0.05 $molid 3.5 30.0 "chain A and resid 1 to 137"   "chain A and resid 1 to 137"
HBpairAnal snr 100.0 90.0 0.05 $molid 3.5 30.0 "chain A and resid 138 to 195" "chain A and resid 138 to 195"

HBresAnal lgn-snr 0 0.0 3.5 30.0 "chain A and resid 1 to 137"   "chain A and resid 138 to 195"
HBresAnal lgn     0 0.0 3.5 30.0 "chain A and resid 1 to 137"   "chain A and resid 1 to 137"
HBresAnal snr     0 0.0 3.5 30.0 "chain A and resid 138 to 195" "chain A and resid 138 to 195"

#DPC
MIanal stBG $molid 100.0 90.0 0.05 4.0 "(resname ARG and name NH1 NH2) or (resname LYS and name NZ) or (resname HSP and name ND1 NE2) or (name N and (same residue as name HT3 H3)) or (resname DPC and name N)" "(resname GLU and name OE2 OE1) or (resname ASP and name OD2 OD1) or name OT1 OT2 O1 O2 or (resname DPC and name O11 O12 O13 O14)"
MIresAnal stBG $molid 0.0 4.0 "(resname ARG and name NH1 NH2) or (resname LYS and name NZ) or (resname HSP and name ND1 NE2) or (name N and (same residue as name HT3 H3)) or (resname DPC and name N)" "(resname GLU and name OE2 OE1) or (resname ASP and name OD2 OD1) or name OT1 OT2 O1 O2 or (resname DPC and name O11 O12 O13 O14)"


MIanal HYDROPHOBIC $molid 40.0 100.0 0.05 4.0 "resname ALA PHE GLY ILE LEU MET PRO VAL TYR and name \"C.*\" and not name C" "resname ALA PHE GLY ILE LEU MET PRO VAL TYR and name \"C.*\" and not name C"
MIresAnal HYDROPHOBIC $molid 0.0 4.0 "resname ALA PHE GLY ILE LEU MET PRO VAL TYR and name \"C.*\" and not name C" "resname ALA PHE GLY ILE LEU MET PRO VAL TYR and name \"C.*\" and not name C"


set molID 0
set refID 1
set seltxt "chain A and name CA"
set frameList {1 2 3}
#set flnmList {rmsfALLapo_crys-2us rmsfALLdpc_crys-2us rmsfALLfar_crys-2us}
set outfile 2usrmsd

set fl [open $outfile w]
set sel [atomselect $molID $seltxt]
set selref [atomselect $refID $seltxt]
foreach uf $frameList {
  $sel frame $uf
  set rmsd($uf) ""
  foreach us [$sel get {x y z}] ur [$selref get {x y z}] {
     set rmsd($uf) "$rmsd($uf) [veclength [vecdist $us $ur]]"
  }
}
foreach uc [$sel get resid] u1 $rmsd(1) u2 $rmsd(2) u3 $rmsd(3) {
  puts $fl "$uc $u1 $u2 $u3"  
}
close $fl



