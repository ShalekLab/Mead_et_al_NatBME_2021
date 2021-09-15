#Analysis_fn.R
#Doug
#Analysis functions used in HTS work
#BY BEN MEAD, DAPHNE JUNE 2017
#VERSION 1.0 (7/19/17)

#-----
#Z_score
z_score = function(tbl, Cols, POSCont){
  tbl.name =  deparse(substitute(tbl))
  Ident.cols = names(tbl)[!names(tbl) %in% Cols]
  scaled = tbl[,Ident.cols]
  
  for(a in Cols){
    if (POSCont == T){
      zscored = tbl[!tbl$Treatment %in% Cont,c(Ident.cols, a)]
      zscored = Controls(tbl,a,POSCont) %>% rbind(zscored)
    }
    if (POSCont ==F){zscored = Controls(tbl,a,POSCont)}
    if (POSCont == "NONE"){zscored = tbl[,c(Ident.cols, a)]}
    values = scale(zscored[,a,drop = F], center = FALSE, scale = apply(zscored[,a,drop = F], 2, sd, na.rm = TRUE))
    zscored[[a]] = values
    names(zscored)[names(zscored) == a] = paste0("Z.", a)
    scaled %<>% left_join(zscored, by = c(Ident.cols))
  }
  cat("\nFound Z scores for", tbl.name, file = pipetxt,append = T)
  return(scaled)
  
}

#Controls: Finds controls for respective populations and tables
#Input: Cont.tbl - table from which to find controls
#       a - assay for which to find controls
#       POSCont - if control is POS (TRUE) or Plate (FALSE)
#Output: Table with controls for indicated assay
Controls = function(Cont.tbl,a,POSCont){
  a = grep(a, names(Cont.tbl),value = T)
  if(length(a) == 0){  cat("\nERROR: Control column not found", file = pipetxt,append = T)}
  #Find control subsets
  if(length(grep("CTG",a))>0){
    Conts = "POS"
    Stims = c(".")}
  if(length(grep("LYZ.S",a))>0){
    Conts = "POS"
    Stims = "NS.S"}
  if(length(grep("LYZ.NS",a))>0){
    Conts = "POS"
    Stims = "^NS."}
  (if (POSCont){
    #find Negative control subset
    POS.tbl = subset(Cont.tbl,Treatment %in% Conts &
                       Stim %in% grep(paste(Stims,collapse = "|"), unique(Cont.tbl$Stim), value = T), 
                     select = c(grep(paste(Ident,collapse = "|"),names(Cont.tbl),value = T),a))}
    
    else{POS.tbl = subset(Cont.tbl, !Treatment %in% c("NC","NONE") & 
                            Stim %in% grep(paste(Stims,collapse = "|"), unique(Cont.tbl$Stim), value = T), 
                          select = c(grep(paste(Ident,collapse = "|"),names(Cont.tbl),value = T),a))})
  if(nrow(POS.tbl) == 0){  cat("\nERROR: Missing control population", file = pipetxt,append = T)}
  return(POS.tbl)
}


#Sum_tbl: Finds mean of primary SSMD analyzed data to combine replicates or match replicate data
#Input: to.sum - data frame with data to summarize
#       Cols - columns to summarize
#       Prim - (TRUE) Primary screen and summarize by indicated column as well as dose and cell type
#            - (FALSE) Replicate screen and summarized by indicated column 
#       By.Col - column by which to summarize (Treatment or Target)
#Output: Data frame with summarized data for indicated columns
Sum_tbl = function(to.sum, Cols, Prim, By.Col){
    tbl.name = deparse(substitute(to.sum))
    Sum.Cols = unique (grep(paste(Cols,collapse="|"), names(to.sum), value=TRUE))
    (if(Prim){groups = c(By.Col,"Dose","Cell.Type")}
    else{groups = By.Col})
    Sum.tbl = to.sum %>% group_by_at(groups)%>% summarize_at(Sum.Cols, funs(mean(.,na.rm = TRUE)))
    cat("\nSummarized",tbl.name,"grouped by",groups, file = pipetxt,append = T)
  return(Sum.tbl)
}

#SSMD_rep_analysis: Finds UMVUE SSMD value (WITH REPLICATES) using Plate as control
#Input:   tbl - table to be analyzed
#         Cols - columns to analyze
#         Cont.tbl - table from which to find controls
#Output: Data frame with Replicate SSMD values
SSMD_rep_UMVUE_analysis = function(tbl, Cols, Cont.tbl){
  tbl = tbl[!tbl$Treatment %in% c("NC", "NONE",Cont),]
  Plates = unique(tbl$Plate)
  SSMD.tbl = tbl[!tbl$Treatment %in% c("NC", "NONE"),c("Treatment"),drop = F] %>% unique()
  for (a in Cols){ #for each assay
    Cont = Controls(Cont.tbl,a,T) 
    Cont$Stim = a 
    Var.tbl = rbind(Cont,tbl[,c(Ident,a)])
    Var.tbl %<>% na.omit() %>% group_by (Treatment, Dose, Cell.Type, Time, Stim) %>% summarise_at(a,funs(mean(.,na.rm = T),var(.,na.rm = T),Rep.n = n())) 
    Var.tbl$mean[Var.tbl$Rep.n <3] = NA
    Var.tbl$var[Var.tbl$Rep.n <3] = NA
    So2 = median(na.omit(Var.tbl$var)) 
    #Calculate SSMD 
    Var.tbl %<>%  mutate(SSMD = (gamma((Rep.n-1)/2)/gamma((Rep.n-2)/2))*sqrt(2/(Rep.n-1))*(mean/sqrt((var+So2)/2)))
    Var.tbl %<>% mutate_if(is.factor, as.character)
    Var.tbl$Rep.n[Var.tbl$Treatment == "POS"] = "Var"
    Var.tbl$Stim[Var.tbl$Treatment == "POS"] = "Var"
    names(Var.tbl)[names(Var.tbl) == "SSMD"] = paste0("repSSMD.",a)
    SSMD.tbl = right_join(SSMD.tbl,Var.tbl[,c("Treatment","Dose", "Cell.Type", "Time", "Stim", "Rep.n",paste0("repSSMD.",a))])
  }
  cat("\nFound Replicate UMVUE SSMD scores with POS as control", file = pipetxt,append = T)
  return(SSMD.tbl)
}

#Save_tbl: Saves CSV of table
#Input:   tbl - table to save
#         S - name to save file as
#Output: NONE - Saves .csv file in Data folder
Save_tbl = function(tbl, S){
  tbl.name = deparse(substitute(tbl))
  name = paste0(Exp.folder,"_",S,".csv")
  write.csv(tbl, paste(Exp.dr,"Analysis_Data",name, sep = "/"),F)
  #append to pipeline
  cat("\nSaved",tbl.name, "as", name, file = pipetxt,append = T)
}

#Crit_FPL_UMVUE: Critical value (SSMD limit) for paired Rep UMVUE SSMD for a FPL of 0.05 with 0.25 cut off as insignificant 
#Input:   tbl - table containing replicate data
#         Col - Define replicates as Biological replicates within a "Treatment" or Treatments within a "Target"
#Output:  Dataframe containing plates and corresponding critical values
#Note - beta2 = 0.25, alpha = 0.05
Crit_FPL_UMVUE= function(tbl,Col){
  #FPL and FNL
  beta2 = 0.25
  tbl %<>% filter(!Treatment %in% c("NC","NONE",Cont))
  #Find v, b, k
  count = as.numeric(unique(tbl[[Col]]))
  (if (any(count<2)) { cat("\nERROR: Replicate count", file = pipetxt,append = T)
    Crit.tbl = data.frame(count, SSMDcrit = NA)}
    else{
      v = count - 1
      b = sqrt(count)
      k = (gamma((count-1)/2)/gamma((count-2)/2))*sqrt(2/(count*(count-1)))
      #set error 
      alpha = 0.05
      #find beta with quantile function
      SSMDcrit = qt((1-alpha),v,ncp = b*beta2)*k
      Crit.tbl = data.frame(count,SSMDcrit,alpha)
    })
  
  x = kable(Crit.tbl, format = "rst", row.names = FALSE)
  cat("\nCritical Values for Replicate screen (",Col,"):\n", file = pipetxt,append = T)
  cat(x,sep = "\n",file = pipetxt,append = T)
  Crit.tbl$alpha = NULL
  return(Crit.tbl)
}
