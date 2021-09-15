#Analysis_Data.R
#Doug
#BEN MEAD, DAPHNE JUNE 2017
#VERSION 1 (11/29/17)

#______
#This script runs data analysis functions
#______
  
  #Find Identifying columns
  Ident = c(names(FC.LOESS.Log.tbl)[!names(FC.LOESS.Log.tbl) %in% Assays])
  
  #Find mean FC for each Treatment and Dose (group biological replicates)
  DoseFC = Sum_tbl(FC.LOESS.Log.tbl, paste0(Assays), Prim = F, c("Treatment","Dose"))
  #Find Zscore of mean FC
  Zscore = z_score(DoseFC, Assays, POSCont = "NONE")
  names(DoseFC)[names(DoseFC) %in% Assays] = paste0("Dose.FC.",names(DoseFC)[names(DoseFC) %in% Assays])
  
  #Find Treatment Rep SSMDs UMVUE (Bio Reps, Keep Treatment and Dose)
  RepSSMD = SSMD_rep_UMVUE_analysis(FC.LOESS.Log.tbl, Assays, FC.LOESS.Log.tbl)
  Stats.run = c("DoseFC", "Zscore","RepSSMD")

  #Combine to final table
  Analyzed.FC.tbl = DoseFC
  for (stat in Stats.run[-1]){
    Analyzed.FC.tbl =  full_join(Analyzed.FC.tbl, get(stat))
  }
  
  #Append Targets
  Analyzed.FC.tbl = right_join(Targ.tbl, Analyzed.FC.tbl)
  
  #Save table
  assign( "Analyzed.FC.tbl",Analyzed.FC.tbl)
  Save_tbl(Analyzed.FC.tbl, "Analyzed.FC.tbl")
    
  rm(RepSSMD, DoseFC,Zscore,stat,Stats.run)

