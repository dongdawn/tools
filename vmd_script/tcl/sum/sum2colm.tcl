#set fileList {dist_A83-B83 dist_A85-B85 dist_A16-B129 dist_A129-B16 dist_chA-chB}
#set fileList {dist_G230-F442 dist_G231-F441 dist_G221-F383 dist_G176-F429 dist_G172-F383 dist_G176-F387 dist_G186-F397 dist_G190-F401 dist_G172-F429}
#set fileList {FGperiL2Angle dist_G231-G148_153 dist_GZres78-G148_153}
#set fileList {dist_waA-sigB dist_waB-sigA}
#set fileList {rmsd2ptNBDz rmsd2ptFGcorez ax_G234@F-8 ax_Gextm3@F-8}
#set fileList {dist_F435-F497 dist_F498-G127 dist_G128-G168 dist_G167-G95 dist_G224-G274 dist_G275-F323 dist_F324-F378 dist_F377-F291 dist_G92-G226 dist_F288-F437}
#set fileList {dist_F205-F252CB dist_F205-F239CB dist_F239-F252CB dist_F92-F205CB dist_F92-F273CB dist_A83-B83CB}
#set fileList {rmsd2ptALL rmsd2ptFG rmsd2ptF rmsd2ptG rmsd2ptFcore rmsd2ptGcore rmsd2ptFGcore rmsd2ptFPL2b rmsd2ptFPL2a rmsd2ptAB rmsd2ptNBD rmsd2ptRD rmsd2ptNBDA rmsd2ptNBDB}
#set fileList {axFtm-axGtm  axGtm-axGich  axFtm-axFich  axFich-axGich axFcore-axGcore axFcored-axGcored axFtm6a-axFtm6b axGtm4a-axGtm4b}
#set fileList {G183@Ftm6_394 G103@Ftm6_394 G107@Ftm6_394 G111@Ftm6_394 F64@Ftm6_394  F67@Ftm6_394  F71@Ftm6_394  F75@Ftm6_394  F79@Ftm6_394 F394@Gtm4_183 F299@Gtm4_183 F303@Gtm4_183 F307@Gtm4_183 G9@Gtm4_183   G13@Gtm4_183  G16@Gtm4_183  G20@Gtm4_183  G24@Gtm4_183  G27@Gtm4_183}
#set fileList {dist_G183-F394 dist_G295-B135 dist_F42-G149 dist_F502-G123 dist_F506-G123 dist_F432-F505 dist_F85-F94 dist_F86-F487 dist_F89-F484 dist_F259-F485 dist_F197-F209 dist_F199-F208 dist_F461-F476 dist_F450-F477 dist_F276-F461 dist_F449-G253 dist_G34-G43 dist_G59-G79 dist_G60-G67 dist_G79-G262 dist_G253-G260 dist_G137-G162 dist_F328-G270 dist_F329-G270 dist_G228-G274 dist_F389-G285 dist_G284-G288 dist_F460-F477}
set fileList {rmsd2inALL rmsd2ptALL rmsd2inNBDA rmsd2inNBDB rmsd2ptNBDA rmsd2ptNBDB rmsd2inNBD rmsd2ptNBD rmsd2inAB rmsd2ptAB rmsd2inFGcore rmsd2ptFGcore rmsd2inGcore rmsd2ptGcore rmsd2inFcore rmsd2ptFcore rmsd2inFG rmsd2ptFG dist_F92-F273CB dist_F92-F205CB dist_F205-F239CB dist_F239-F252CB dist_F205-F252CB dist_G183-F394 dist_G183-F394 dist_G176-F429 dist_G221-F383 dist_G230-F442 dist_G231-F441 dist_NBDA-NBDB dist_waA-sigB dist_sigA-waB}
set docuList {"W:/weng/inE2ptE/1"      inE2ptE  1    0 0 \
              "W:/weng/inE2ptE/2"      inE2ptE  2    0 0 \
              "W:/weng/inE2ptE/3"      inE2ptE  3    0 0 \
              "W:/weng/inE2ptE/4"      inE2ptE  4    0 0 \
              "W:/weng/inE2ptE/5"      inE2ptE  5    0 0 \
              "W:/weng/inE2ptE/5ns1"   inE2ptE  5ns1 0 0 \
              "W:/weng/inE2ptE/5ns2"   inE2ptE  5ns2 0 0 \
              "W:/weng/ptE2inE/1"      ptE2inE  1    0 0 \
              "W:/weng/ptE2inE/2"      ptE2inE  2    0 0 \
              "W:/weng/ptE2inE/3"      ptE2inE  3    0 0 \
              "W:/weng/ptE2inE/4"      ptE2inE  4    0 0 \
              "W:/weng/ptE2inE/5"      ptE2inE  5    0 0 \
              "W:/weng/ptE2inE/5ns1"   ptE2inE  5ns1 0 0 \
              "W:/weng/ptE2inE/5ns2"   ptE2inE  5ns2 0 0 \
             }

foreach uf $fileList {
 # read file #
  set hbList ""
  set idx 0
  foreach {ud u45 uro ude uwm} $docuList {
    set labl($idx) "$u45 $uro $ude $uwm"
    if {[file exists ${ud}/${uf}]} then {
    set rfl [open ${ud}/${uf} r]
    while {[gets $rfl line]>3} {
      lappend hbList [lindex $line 0]
      set occ("[lindex $line 0]_${idx}") [lindex $line 1]
    }
    close $rfl
    }
    incr idx
  }
  set ndx $idx
 # output #
  set wfl [open ${uf}SUM w]
  set out45 "-"
  set outro "-"
  set outde "-"
  set outwm "-"
  for {set idx 0} {$idx < $ndx} {incr idx} {
    set out45 "$out45 [lindex $labl($idx) 0]"
    set outro "$outro [lindex $labl($idx) 1]"
    set outde "$outde [lindex $labl($idx) 2]"
    set outwm "$outwm [lindex $labl($idx) 3]"
  }
  puts $wfl $out45
  puts $wfl $outro
  puts $wfl $outde
  puts $wfl $outwm
  foreach uhb [lsort -unique -integer $hbList] {
    set outtxt $uhb
    for {set idx 0} {$idx < $ndx} {incr idx} {
      if {[info exists occ(\"${uhb}_${idx}\")]} then {
        set outtxt "$outtxt $occ(\"${uhb}_${idx}\")"
      } else {
        set outtxt "$outtxt 0"
      }
    }
    puts $wfl $outtxt
  }
  close $wfl
  unset occ
}


#set docuList {"E:/weng/inMalFGK1"         MalFGK  apo1    0 0  \
#              "E:/weng/inMalFGK2"         MalFGK  apo2    0 0  \
#              "E:/weng/inMalFGK3"         MalFGK  apo3    0 0  \
#              "E:/weng/inMalFGKatp1"      MalFGK  atp1    0 0  \
#              "E:/weng/inMalFGKatp2"      MalFGK  atp2    0 0  \
#              "E:/weng/inMalFGKatp3"      MalFGK  atp3    0 0  \
#              "E:/weng/inMalFGK-MalE1"    MalFGKE apo1    0 0  \
#              "E:/weng/inMalFGK-MalE2"    MalFGKE apo2    0 0  \
#              "E:/weng/inMalFGK-MalE3"    MalFGKE apo3    0 0  \
#              "E:/weng/inMalFGK-MalEatp1" MalFGKE atp1    0 0  \
#              "E:/weng/inMalFGK-MalEatp2" MalFGKE atp2    0 0  \
#              "E:/weng/inMalFGK-MalEatp3" MalFGKE atp3    0 0  \
#             }

#set docuList {"W:/weng/in2pt/1"  in2pt  1  0 0 \
#              "W:/weng/in2pt/2"  in2pt  2  0 0 \
#              "W:/weng/in2pt/3"  in2pt  3  0 0 \
#              "W:/weng/in2pt/4"  in2pt  4  0 0 \
#              "W:/weng/in2pt/5"  in2pt  5  0 0 \
#              "W:/weng/pt2in/1"  pt2in  1  0 0 \
#              "W:/weng/pt2in/2"  pt2in  2  0 0 \
#              "W:/weng/pt2in/3"  pt2in  3  0 0 \
#              "W:/weng/pt2in/4"  pt2in  4  0 0 \
#              "W:/weng/pt2in/5"  pt2in  5  0 0 \
#             }

