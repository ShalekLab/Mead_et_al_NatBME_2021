#Pipeline.R
#PIPELINE RUN
#Optimized 1w
#BEN MEAD, DAPHNE SZE JUNE 2017
#VERSION 2.0 (6/25/17)

#______
#This script will call all other scripts & fns to assemble a complete pipeline of quality control tests
#______

#Install missing packages
packages = c("readxl", "dplyr","tidyr", "ggplot2", "ggrepel","RColorBrewer","gridExtra","magrittr","knitr","reshape2","fdrtool")
to.install = packages[!(packages %in% installed.packages()[,"Package"])]
if(length(to.install)) install.packages(to.install)
#Load Libraries
lapply(packages, require, character.only = TRUE)
rm(packages, to.install)

#Get functions, experiment, and pipeline.txt directory
fcn.dr = getwd()
Exp.dr = paste(fcn.dr,Exp.folder, sep = "/")
pipetxt = paste(Exp.dr,paste0(Exp.folder,"_Quality_pipe.txt"),sep = "/")


#Load experiment details and global variables 
Exp.det = read_excel(paste(Exp.dr,"RAW_Data/Experiment_Details.xlsx",sep="/"),col_names=T)
Assays = na.omit(unlist(Exp.det$Assays))
Cont  = na.omit(unlist(Exp.det$Cont))
Dates = unique(na.omit(unlist(Exp.det$Dates)))
C.Std = na.omit(unlist(Exp.det$C.Std))[1]
Target.file = na.omit(unlist(Exp.det$Target.file))[1]

#Initalize pipeline.txt
if (!file.exists(pipetxt)){
  cat("", file = pipetxt)
}
if(length(grep("Assays performed:",readLines(pipetxt))) == 0){
  #save to pipeline
  cat(Exp.folder, "pipeline",file = pipetxt)
  cat("\nTimes pipeline was run:\n", 0 , file = pipetxt, append = T)
  cat("\nAssays performed: ",Assays , file = pipetxt, append = T)
  cat("\nControls used: ",Cont , file = pipetxt,append = T)
  cat("\nDates of Experiments: ",Dates , file = pipetxt,append = T)
}
#save run iteration
pipe.ran = as.integer(readLines(pipetxt)[3])
cat("\n\nPipeline run:",pipe.ran +1, "\ton:", as.character.Date(Sys.time()),file = pipetxt,append = T)


#Import Target lists 
Targ.tbl = read_excel(paste(Exp.dr, "RAW_Data",Target.file,sep="/"),col_names=T, sheet = "Old") 
new.Targ = read_excel(paste(Exp.dr, "RAW_Data",Target.file,sep="/"),col_names=T, sheet = "New")
#Combine and update old and new Target lists
Targ.tbl %<>% left_join(new.Targ, by = "Treatment")
changes = Targ.tbl %>% filter(Target.x != Target.y)
changes.tbl = Targ.tbl %>% filter(Target.x %in% changes$Target.x| Target.y %in% changes$Target.y)
Targ.tbl$Target.x[!is.na(Targ.tbl$Target.y) & !is.na(Targ.tbl$Target.x) &!Targ.tbl$Target.x == Targ.tbl$Target.y] =
  Targ.tbl$Target.y[!is.na(Targ.tbl$Target.y) & !is.na(Targ.tbl$Target.x) &!Targ.tbl$Target.x == Targ.tbl$Target.y]
#Finalize Target table
Targ.tbl$Treatment = factor(Targ.tbl$Treatment)
Targ.tbl %<>% mutate(Target = paste(Target.x, Activity)) %>% select(Treatment, Target, Pathway)
cat("\nImported Target list", file = pipetxt,append = T)

rm(Target.file)
#_____

#Plots functions
source("./R_scripts_fns/Manu_Plots.R")

#Import and convert raw data to tables
#IMPORT from .xlsx files to Clist & assay arrays
source("./R_scripts_fns/Import.R")

#CONVERT raw arrays to B-scored and NC subtracted final tables
source("./R_scripts_fns/Convert.R")

#_____
#Load analysis functions
source("./R_scripts_fns/Analysis_fn.R")
#Select tables to run analysis
cat("\n\nPerforming Data Analysis for: FC.Full.tbl" , file = pipetxt,append = T)
#Run analysis on Data
source("./R_scripts_fns/Analysis_Data.R")

#Select tables to find hits 
cat("\n\nPerforming Replicate Hit Analysis for: FC.Full.tbl", file = pipetxt,append = T)
#Plot treatment, target and dose response graphs
source("./R_scripts_fns/Analysis_Hits_Rep.R")

#Plot Plots

FPL.Range()
RepSSMD_plot()
Zscore_plot()
Full_Dose()
Time_Dose()
if(LOESS.Norm == T){LOESS_plates()}
dual_flashlight()

#update pipeline.txt
pipe.ran = pipe.ran+ 1
pipeline = readLines(pipetxt,-1)
pipeline[3] = as.character(pipe.ran)
writeLines(pipeline,pipetxt,sep = "\n")

#append pipe to README
# pipeline = readLines("pipeline.txt",-1)
# writeLines(pipeline,"README.txt",sep = "\n")

#clear global environment
# rm(list=ls())
