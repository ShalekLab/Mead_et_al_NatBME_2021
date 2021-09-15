#Analysis_Hits_Rep.R

  #Select tables
  Full.tbl = Analyzed.FC.tbl

  #Columns for Replicate SSMD
  Rep.cols = grep("_Rep",names(Full.tbl),value = T)
  
  #find critical values for RepSSMD
  RepCrit = Crit_X_UMVUE(Full.tbl,"Rep.n")
  RepCrit %<>% rename(RepCrit = SSMDcrit) %>% mutate(Rep.n = as.character(count))
  RepCrit$count = NULL
  Full.tbl %<>% full_join(RepCrit)
  
  #Filter  based on Rep SSMD cutoff
  Filt.tbl = Full.tbl %>% filter(SSMD.LYZ.S_Rep > RepCrit  & SSMD.LYZ.NS_Rep > RepCrit) 
  Filt.tbl$RepCrit = NULL
  
  #Append Targets
  Filt.tbl = right_join(Targ.tbl, Filt.tbl)
  
  assign( "Rep_Hits.tbl",Filt.tbl)
  Save_tbl(Filt.tbl, paste0("Rep_hits"))
  
  #find critical value for Z score to 10%
  Zcrit = qnorm(0.9)
  cat("\nCritical Values for Z scored FC (10%):", Zcrit ,file = pipetxt,append = T)
  
  #Filter  based on Zscore cutoff
  Filt.tbl %<>% filter(Z.LYZ.S > Zcrit & Z.LYZ.NS > Zcrit)
  assign("Z_Rep_Hits.tbl",Filt.tbl)
  Save_tbl(Filt.tbl, paste0("Z_Rep_hits"))
  
  
  rm(Rep.cols,Zcrit,RepCrit,Full.tbl,Filt.tbl)
  