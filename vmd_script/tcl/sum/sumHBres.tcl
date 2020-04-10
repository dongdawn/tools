#set fileList {chEa-FGHISpairP chEb-FGHISpairP}
#set docuList {"E:/weng/inMalFGK1"         MalFGK  apo1    0 0 \
#              "E:/weng/inMalFGK2"         MalFGK  apo2    0 0 \
#              "E:/weng/inMalFGK3"         MalFGK  apo3    0 0 \
#              "E:/weng/inMalFGKatp1"      MalFGK  atp1    0 0 \
#              "E:/weng/inMalFGKatp2"      MalFGK  atp2    0 0 \
#              "E:/weng/inMalFGKatp3"      MalFGK  atp3    0 0 \
#              "E:/weng/inMalFGK-MalE1"    MalFGKE apo1    0 0 \
#              "E:/weng/inMalFGK-MalE2"    MalFGKE apo2    0 0 \
#              "E:/weng/inMalFGK-MalE3"    MalFGKE apo3    0 0 \
#              "E:/weng/inMalFGK-MalEatp1" MalFGKE atp1    0 0 \
#              "E:/weng/inMalFGK-MalEatp2" MalFGKE atp2    0 0 \
#              "E:/weng/inMalFGK-MalEatp3" MalFGKE atp3    0 0 \
#             }

#set fileList {C2_275-308RESP C2_276-308RESP C2_279-308RESP C2_278-306RESP C2_252-303RESP C2_252-305RESP C2_252-306RESP C2_253-279RESP C2_205-302RESP C2_205-305RESP}
#set docuList {"E:/yanxin/research/ABC_data/MetNI/my_work/6.md_with_3tuz/analysis"  holo met 0 0 \
#              "E:/yanxin/research/ABC_data/MetNI/my_work/4.md_with_MetNI/MetN"     apo   -  0 0 \
#              "E:/yanxin/research/ABC_data/MetNI/my_work/4.md_with_MetNI/MetN+Met" apo  met 0 0 \
#              "E:/yanxin/research/ABC_data/MetNI/my_work/6.md_with_3tuz/MetN/new"  holo  -  0 0 \
#             }
             
set fileList {stBGHISpairP}
set docuList {"L:/ykt6/ykt6apo"                           apo - 0 0 \
              "L:/ykt6/ykt6dpc（输出误写成apo，应为dpc）" dpc - 0 0 \
              "L:/ykt6/ykt6far"                           far - 0 0 \
             }
set spc " - "

foreach uf $fileList {
 # read HBfile #
  set hbList ""
  set idx 0
  foreach {ud u45 uro ude uwm} $docuList {
    set labl($idx) "$u45 $uro $ude $uwm"
    if {[file exists ${ud}/${uf}]} then {
      set rfl [open ${ud}/${uf} r]
      while {[gets $rfl line]>10} {
        lappend hbList "[lrange $line 0 2] - [lrange $line 4 6]"
        set occ("[lrange $line 0 2]$spc[lrange $line 4 6]_${idx}") [lindex $line 7]
      }
      close $rfl
    } else {
      puts "${ud}/${uf} NOT exists!"  
    }
    incr idx
  }
  set ndx $idx
 # output #
  set wfl [open ${uf}SUM w]
  set out45 "- - - - - - -"
  set outro "- - - - - - -"
  set outde "- - - - - - -"
  set outwm "- - - - - - -"
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
  foreach uhb [lsort -unique $hbList] {
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
  if {[info exists occ]} {unset occ}
}

