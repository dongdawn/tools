proc calcrmsf {selcalcList selovaList flnmList refID molID initF lastF} {
  foreach usc $selcalcList uso $selovaList flnm $flnmList {
    set fl [open ${flnm}_${initF}-${lastF} w]
    alignTRAJ $refID $molID $uso
    set sel [atomselect $molID $usc]
    foreach unt [$sel get index] chres [$sel get {chain resid}] {
      set tmp [measure rmsf [atomselect $molID "index $unt"] first $initF last $lastF]
      puts $fl "$chres $tmp"
    }
    close $fl
  }
}
