#Convert.R
#BEN MEAD, DAPHNE SZE JUNE 2017
#VERSION 2.0 (6/22/17)
#_______
#This script converts the RAW plate arrays into a B-scored table, raw data corrected for NC controls, LOESS normalized table
#_______
#load fns
source("./R_scripts_fns/Convert_fn.R")
Experiments = vector("list",length(Dates))

#Compute raw values for all plates
for (i in 1:length(Experiments)){
  Experiments[[i]] = Convert(as.character(Dates[i]))
}
#Combine EXPTS
RAW.Full.tbl=do.call(bind_rows,Experiments)
RAW.Full.tbl %<>% mutate_at(c("Treatment","Dose","Cell.Type","Stim","Time"),funs(factor(.)))
RAW.Full.tbl %<>% filter(!is.na(RAW.Full.tbl$CTG) & !Treatment =="REMOVE")

Save_tbl(RAW.Full.tbl,"Raw_data")
#Append to pipeline.txt
cat("\n\t RAW.Full.tbl", file = pipetxt,append = T)

#Find Log data
Log.Raw.tbl = RAW.Full.tbl
for(a in Assays){
  Log.Raw.tbl[[a]] = log10(RAW.Full.tbl[[a]])
}
cat("\nTransformed raw data to Log10 space", file = pipetxt ,append = T)
Save_tbl(Log.Raw.tbl,"Log_Norm_data")

Ident = c(names(Log.Raw.tbl)[!names(Log.Raw.tbl) %in% Assays])


(if(LOESS.Norm == T){
  run = "LOESS"
  #Find Normalized LOESS data
  LOESS.Log.tbl = LOESS_norm(Log.Raw.tbl)
  Save_tbl(LOESS.Log.tbl,"LOESS_Norm_data")
  #Find FC data
  FC.Full.tbl = SubLog_FC(LOESS.Log.tbl,Assays,T,LOESS.Log.tbl)
  }
  else{
  run = "FC"
  #Find FC data
  FC.Full.tbl = SubLog_FC(Log.Raw.tbl,Assays,T,Log.Raw.tbl)
  
  })
  Save_tbl(FC.Full.tbl, "Plate_Norm_data")

  Distributions(run)
  Bio.fn(run)
  Bio.plate.fcn(run)
  POS_dist(run)
  
# rm(Experiments,i,Clists,Alists,Qlists,Dates,a, Ident)
#rm(RAW.Full.tbl,Log.Raw.tbl,LOESS.Log.tbl)

