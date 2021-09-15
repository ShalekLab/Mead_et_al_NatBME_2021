#Analysis_Data.R
#Doug
#BEN MEAD, DAPHNE JUNE 2017
#VERSION 1 (11/29/17)

#______
#This script runs data analysis functions
#______
for (tbl in tbls){
  #Table to analyze
  m.FC.LOESS.Log.tbl = Data.list[[tbl]]
  cat("\n\nData Analysis for: ", tbl , file = pipetxt, append = T)
  
  #Find Identifying columns
  Ident = c(names(m.FC.LOESS.Log.tbl)[!names(m.FC.LOESS.Log.tbl) %in% Assays])
  
  DoseFC = Sum_tbl(m.FC.LOESS.Log.tbl, paste0(Assays), Prim = F, c("Treatment","Dose"))
  Zscore = z_score(DoseFC, Assays, POSCont = "NONE")
  names(DoseFC)[names(DoseFC) %in% Assays] = paste0("Dose.FC.",names(DoseFC)[names(DoseFC) %in% Assays])
  
  #Find Treatment Rep SSMDs UMVUE (Bio Reps, Keep Treatment and Dose)
  PrimSSMD = SSMD_primary_analysis(m.FC.LOESS.Log.tbl, Assays, POSCont = F, m.FC.LOESS.Log.tbl)
  Stats.run = c("DoseFC", "Zscore","PrimSSMD")
  
  #Combine to final table
  Full.tbl = DoseFC
  for (stat in Stats.run[-1]){
    Full.tbl =  full_join(Full.tbl, get(stat))
  }
  
  #Append Targets
  Full.tbl = right_join(Targ.tbl, Full.tbl)
  
  assign( paste0("Full_",tbl),Full.tbl)
  Save_tbl(Full.tbl,paste0("Full_",tbl,".tbl"))
  
  rm( PrimSSMD, Zscore, DoseFC, Full.tbl, stat, Stats.run)
}
rm(tbl,tbls)
