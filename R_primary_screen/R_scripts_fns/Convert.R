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
Raw.tbl=do.call(bind_rows,Experiments)
Raw.tbl %<>% mutate_at(c("Treatment","Dose","Cell.Type","Stim"),funs(factor(.)))
#Append to pipeline.txt
cat("\nConverted to raw data to table with columns:\n",names(Raw.tbl), file = pipetxt ,append = T)
Save_tbl(Raw.tbl,"Raw_data")

#Find Log data
Log.Raw.tbl = Raw.tbl
for(a in Assays){
  Log.Raw.tbl[[a]] = log10(Raw.tbl[[a]])
}
cat("\nTransformed raw data to Log10 space", file = pipetxt ,append = T)
Save_tbl(Log.Raw.tbl,"Log_Norm_data")
#Find Normalized LOESS data
LOESS.Log.tbl = LOESS_norm(Log.Raw.tbl)
Save_tbl(LOESS.Log.tbl,"LOESS_Norm_data")

Ident = c(names(LOESS.Log.tbl)[!names(LOESS.Log.tbl) %in% Assays])

#Find FC data
FC.LOESS.Log.tbl = SubLog_FC(LOESS.Log.tbl,Assays,F,LOESS.Log.tbl)
Save_tbl(FC.LOESS.Log.tbl, "Plate_Norm_data")

rm(Experiments,i,Clists,Alists, Plists,Qlists,File.Type,Dates,a, Ident)
rm(Raw.tbl, Log.Raw.tbl)

rm(list = lsf.str())
