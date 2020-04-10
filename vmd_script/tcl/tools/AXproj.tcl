proc AXproj {molID selaxtxt selotxt outfileList selpntList} {
  set selax [atomselect $molID $selaxtxt]
  set selo [atomselect $molID $selotxt]
  set fls ""
  set selpnts ""
  foreach uf $outfileList us $selpntList {
    set usel [atomselect $molID $us]
    if {[$usel num] >0} {
      lappend selpnts $usel
      lappend fls [open $uf w]
    } else {
      puts "NO $us !"
    }
  }
  set nframe [molinfo $molID get numframes]
  for {set iframe 0} {$iframe < $nframe} {incr iframe} {
    $selax frame $iframe
    $selo frame $iframe
    set ax [lindex [lindex [measure inertia $selax] 1] 2]
    set oxyz [measure center $selo weight mass]
    foreach uf $fls us $selpnts {
      $us frame $iframe
      set uxyz [measure center $us weight mass]
      set tp [ vecdist $uxyz [vecadd $oxyz [vecscale [vecdot [vecsub $uxyz $oxyz] $ax] $ax]] ]
      puts $uf "$iframe $tp"
    }
  }
  foreach uf $fls {close $uf}
}
#    molID                selaxtxt                               selotxt                    outfileList   selpntList
#AXproj 0 "chain F and resid 381 to 394 and name CA" "chain F and resid 394 and name CA"  \
       { G183@Ftm6_394                          G103@Ftm6_394           \
         G107@Ftm6_394                          G111@Ftm6_394           \
         F64@Ftm6_394                           F67@Ftm6_394            \
         F71@Ftm6_394                           F75@Ftm6_394            \
         F79@Ftm6_394                     } \
       {"chain G and resid 183 and name CA"     "chain G and resid 103 and name CA"  \
        "chain G and resid 107 and name CA"     "chain G and resid 111 and name CA"  \
        "chain F and resid  64 and name CA"     "chain F and resid  67 and name CA"  \
        "chain F and resid  71 and name CA"     "chain F and resid  75 and name CA"  \
        "chain F and resid  79 and name CA" }

#AXproj 0 "chain G and resid 172 to 183 and name CA" "chain G and resid 183 and name CA" \
       { F394@Gtm4_183                          F299@Gtm4_183                \
       	 F303@Gtm4_183                          F307@Gtm4_183                \
         G9@Gtm4_183                            G13@Gtm4_183                 \
         G16@Gtm4_183                           G20@Gtm4_183                 \
         G24@Gtm4_183                           G27@Gtm4_183              }  \
       {"chain F and resid 394 and name CA"     "chain F and resid 299 and name CA"  \
       	"chain F and resid 303 and name CA"     "chain F and resid 307 and name CA"  \
        "chain G and resid   9 and name CA"     "chain G and resid  13 and name CA"  \
        "chain G and resid  16 and name CA"     "chain G and resid  20 and name CA"  \
        "chain G and resid  24 and name CA"     "chain G and resid  27 and name CA" }


