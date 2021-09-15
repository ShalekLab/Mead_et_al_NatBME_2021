#Analysis_Data.R
#Doug
#BEN MEAD, DAPHNE JUNE 2017
#VERSION 1 (11/29/17)

#______
#This script runs data analysis functions
#______
  #Remove No Cell, NONE and Remove Treatments  
  FC.Full.tbl %<>% filter(!Treatment %in% c("NC","NONE","REMOVE") & Time %in% c("FULL","NONE"))

  #Find Identifying columns
  Ident = c(names(FC.Full.tbl)[!names(FC.Full.tbl) %in% Assays])
  
  #Find mean FC for each Treatment and Dose (group biological replicates)
  DoseFC = Sum_tbl(FC.Full.tbl, paste0(Assays), Prim = F, c("Time","Treatment","Dose"))
  #Find Zscore of mean FC
  Zscore = z_score(DoseFC, Assays, POSCont = "NONE")
  names(DoseFC)[names(DoseFC) %in% Assays] = paste0("Dose.FC.",names(DoseFC)[names(DoseFC) %in% Assays])
  
  #Find Treatment Rep SSMDs UMVUE (Bio Reps, Keep Treatment and Dose)
  RepSSMD = SSMD_rep_UMVUE_analysis(FC.Full.tbl, Assays, FC.Full.tbl)
  Stats.run = c("DoseFC", "Zscore","RepSSMD")

  Analyzed.FC.tbl = DoseFC
  for (stat in Stats.run[-1]){
    Analyzed.FC.tbl =  full_join(Analyzed.FC.tbl, get(stat))
  }

  assign( "Analyzed.FC.tbl", Analyzed.FC.tbl)
  Save_tbl(Analyzed.FC.tbl, "Analyzed.FC.tbl")
  
  #Select relevant Data for Dose plots
  Dose.data = FC.Full.tbl %>% filter(Time == "FULL" & !Treatment %in% c("Cont","NC")) %>% select("Treatment","Dose","Plate","Quad", Assays)  %>% rename(Rep = Plate) %>% mutate(Dose = as.numeric(as.character(Dose)))
  #Find Doses per treatment
  Doses.val = distinct(Dose.data[,c("Treatment","Dose")]) %>% group_by(Treatment) %>% filter(Dose == max(Dose)) %>% mutate(x2 = Dose/2, x4 = Dose/4, x8 = Dose/8)
  Save_tbl(Doses.val,"Treatment.Doses")
  #Set Dose to generic 
  Doses.val = melt(Doses.val, id.vars = "Treatment") %>% rename(Dose = value)
  Dose.data %<>% left_join(Doses.val) %>% select(-Dose) %>% rename(Dose = variable)

  #Remove Date identifier
  Dose.data$Rep = sub(".* ","",Dose.data$Rep)
  #Formate table
  for (a in Assays){
    Dose.tbl = dcast(Dose.data[,c("Treatment","Rep","Quad","Dose",a)],Treatment + Rep + Quad ~ Dose, value.var = a) 
    Save_tbl(Dose.tbl, paste0(a,"_Dose_Data"))
  }
    
  rm(RepSSMD, DoseFC,stat,Stats.run)

