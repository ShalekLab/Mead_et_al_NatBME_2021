#Convert_fn.R
#BEN MEAD, DAPHNE SZE June 2017
#VERSION 2 (6/22/17)


#_____
#Convert:calls functions to convert arrays into a table
#Input: 
#     date: date of experiment (used for calling files - MD) 
#     BS: if data should be B scored
#     NC: if data should be NC subtracted
#Output: Merged table

Convert=function(date){
  Qlist = Qlists[[date]]
  #Build NC subtracted arrays
  #read Compound, Dose and POS stim lists
  Clist = Clists[[date]]
  Dlist = Dlists[[date]]
  POSlist = POSlists[[date]]
  #Data lists
  Alist = Alists[[date]]
  #Split into quadrants
  Q.B.list =  vector("list",length(Assays))
  names(Q.B.list) = Assays
  for (a in Assays){
    Q.B.list[[a]] =  Arr_to_Q(Alist[[a]],Qlist)
  }
  rm(a)
  #convert data arrays to single table
  B.tbls =  vector("list",length(Assays))
  names(B.tbls) = Assays
  for (a in Assays){
    B.tbls[[a]]= Arr_to_tbl(Q.B.list[[a]],Clist,a,date, Qlist, POSlist, Dlist)
  }
  #Merge table (check if 2 or 3 assays performed)
  Merged.tbl=Mrg_tbl(B.tbls)
  rm(B.tbls,Clist,Q.B.list)
  #return merged table
  return(Merged.tbl)
}

#______
#Arr_to_Q: Breaks 16x24 (384) well plates into 8x12 (96) sets based on a 'quadrant' label
#Input: B.arr: a 384w array containing assay data
#       Qlist: a list of vectors containing quadrant info per below:
#               (1)(2)
#               (3)(4)
#Output: Q.B.arr: list of quadrant split arrays for each  384w plate
Arr_to_Q=function(B.Arr,Qlist) {
  #dimentions and names for 384 well plate
  d=dim(B.Arr)
  dn=dimnames(B.Arr)
  #Dimentions for quadrant (96w plate)
  dQ=c(d[1]/2,d[2]/2,4)
  #initalize list of plates
  Q.B.Arr = vector("list", d[3])
  names(Q.B.Arr) = dn[[3]]
  for(i in 1:d[3]) {      # 384w plate
    #initialize array for quadrants
    Q.B.Arr[[i]]=array(c(0),dQ,dimnames = list(NULL,NULL,Qlist[[i]]))
    for(j in 1:dQ[3]) {    # quadrant
      for(k in 1:dQ[2]) {  # column
        for(l in 1:dQ[1]){ # row
          if (j==1) {
            Q.B.Arr[[i]][l,k,j]=B.Arr[2*l-1,2*k-1,i]
          } else if (j==2) {
            Q.B.Arr[[i]][l,k,j]=B.Arr[2*l-1,2*k,i]
          } else if (j==3) {
            Q.B.Arr[[i]][l,k,j]=B.Arr[2*l,2*k-1,i]
          } else
            Q.B.Arr[[i]][l,k,j]=B.Arr[2*l,2*k,i]
        }
      }
    }
  }
  return(Q.B.Arr) 
}
#______
#Arr_to_tbl: change array into formatted data table 
#Input: 
#     Q.B.Arr: list of quadrant split 8x12 arrays (96 well) for each 384w plate
#     Clist: 8x12 (96) well compound/treatment structure for screen
#     Asssay: name of Assay
#     Date: Date of experiment
#     Qlist: list of vectors containing quadrant data for each 384w plate
#Output: formatted table
Arr_to_tbl=function(Q.B.Arr,Clist,Assay,date,Qlist, POSlist, Dlist) {
  #Convert compound list to formatted table
  #name columns with A to L
  dimnames(Clist)[2]=list(LETTERS[1:12])
  #Create Identifying columns
  Clist= Clist %>%
    as.data.frame.table(stringsAsFactors = F) %>%
    mutate(Index=paste(Var1,Var2)) %>%
    mutate(ID=paste(Var1,Var2,Var3))
  #Rename column
  names(Clist)[names(Clist) == "Freq"] ="Treatment"
  
  Dlist = Dlist %>%
    as.data.frame.table(stringsAsFactors = F) %>%
    mutate(Index=paste(Var1,Var2)) %>%
    mutate(ID=paste(Var1,Var2,Var3))
  #Rename column
  names(Dlist)[names(Dlist) == "Freq"] ="Dose"
  
  POSlist = POSlist %>%
    as.data.frame.table(stringsAsFactors = F)%>%
    mutate(Index=paste(Var1,Var2)) %>%
    mutate(ID=paste(Var1,Var2,Var3))
  #Rename column
  names(POSlist)[names(POSlist) == "Freq"] ="POS.Stim"
  
  #Join tables
  Clist %<>% full_join(Dlist, by = c("Var1", "Var2", "Var3", "Index", "ID")) %>% full_join(POSlist,by = c("Var1", "Var2", "Var3", "Index", "ID"))
  
  #Select relevant colums
  Clist %<>% select(c("Treatment","Dose","POS.Stim","ID","Index"))
  
  #data array to formatted table
  B.list = vector("list",length(Q.B.Arr))
  for (i in 1:length(Q.B.Arr)){
    B.list[[i]] = as.data.frame.table(Q.B.Arr[[i]],stringsAsFactors = F)
    names(B.list[[i]]) = c("row","col","Quad",Assay)
    B.list[[i]]$Plate = names(Q.B.Arr)[i]
  }
  B.tbl=do.call(bind_rows,B.list)
  
  #Find ID Var
  B.tbl = mutate(B.tbl,ID=paste(row,col,Plate))
  #Add row and columns for 384 well plate
  #Convert rows and cols into numbers
  B.tbl$row = match(B.tbl$row,LETTERS[1:8])
  B.tbl$col = match(B.tbl$col,LETTERS[1:12])
  for (p in names(Qlist)){ #for each plate
    #find name of Quad in each each row and col
    row1 = c(Qlist[[p]][c(1,2)])
    row2 = c(Qlist[[p]][c(3,4)])
    col1 = c(Qlist[[p]][c(1,3)])
    col2 = c(Qlist[[p]][c(2,4)])
    #Find 384 row and col number
    B.tbl$row[B.tbl$Plate == p & B.tbl$Quad %in% row1] = (2*B.tbl$row[B.tbl$Plate == p & B.tbl$Quad %in% row1])-1
    B.tbl$row[B.tbl$Plate == p & B.tbl$Quad %in% row2] = 2*B.tbl$row[B.tbl$Plate == p & B.tbl$Quad %in% row2]
    B.tbl$col[B.tbl$Plate == p & B.tbl$Quad %in% col1] = (2*B.tbl$col[B.tbl$Plate == p & B.tbl$Quad %in% col1])-1
    B.tbl$col[B.tbl$Plate == p & B.tbl$Quad %in% col2] = 2*B.tbl$col[B.tbl$Plate == p & B.tbl$Quad %in% col2]
  }
  
  #rename & format tabl3
  #Changes made here
  B.tbl = left_join(Clist,B.tbl,by = "ID")
  B.tbl %<>% mutate_if(is.factor, as.character) 
  B.tbl %<>% filter(!Treatment == "NONE")
  B.tbl$ID = NULL
  B.tbl = separate(B.tbl, Quad, c("Quad","Time","Cell.Type","Stim"), "[ ]")
  B.tbl=mutate(B.tbl,Plate=paste(date,Plate))
  #Fix Controls
  B.tbl$Stim[B.tbl$Stim == "POS"] = as.character(B.tbl$POS.Stim[B.tbl$Stim == "POS"])
  B.tbl$Treatment[B.tbl$Time == "NONE" & B.tbl$Treatment != "NC"] = "POS"
  B.tbl$Dose[B.tbl$Time == "NONE" & B.tbl$Treatment != "NC"] = 0
  B.tbl$POS.Stim = NULL
  
  B.tbl$Stim[B.tbl$Treatment == "NC"] = "NS.NS"
  
  #return tables
  return(B.tbl)
}

#______
#Mrg_tbl:combine measurements into single data table
#Input: List of tables to merge
#Output: merged table with duplicate rows removed
Mrg_tbl=function(Arr.list) {
  #Join tables
  mtbl = Arr.list[[1]]
  for (a in Assays[-1]){
    mtbl = left_join(mtbl,Arr.list[[a]],by = c("Treatment", "Index","Quad", "row", "col", "Dose", "Cell.Type", "Stim", "Plate", "Time"))
  }
  #remove empty Control rows
  mtbl=mtbl[c(names(mtbl)[!names(mtbl) %in% Assays],Assays)]
  return(mtbl)
}

#SubLog_FC: Finds the fold change over the control per Plate (Val - mean(Control))
#Input: tbl - table to be analyzed
#       Cols - columns to analyze
#       POSCont - if control is POS (TRUE) or Plate (FALSE)
#       Cont.tbl - table from which to find controls
#Output: data frame with POS FC columns
SubLog_FC = function(tbl, Cols, POSCont, Cont.tbl){
  p.tbl = tbl[!tbl$Treatment %in% c("NONE"),]
  for (c in Cols){ # for each assay
    #Find control populations
    POS.tbl = Controls(Cont.tbl,c,POSCont)
    for (p in unique(POS.tbl$Plate)){ # for each plate
      #find the mean of the control for the plate
      POS = median(POS.tbl[[c]][POS.tbl$Plate == p],na.rm = T)
      #Calculate  FC
      p.tbl[[c]][p.tbl$Plate == p] = p.tbl[[c]][p.tbl$Plate == p] - POS
    }
  }
  #append to pipeline.txt
  cat("\nFound fold change of log values relative to median of",ifelse(POSCont,"POS","plate"), file = pipetxt,append = T)
  return(p.tbl)
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

#Save_tbl: Saves CSV of table
#Input:   tbl - table to save
#         S - name to save file as
#Output: NONE - Saves .csv file in Data folder
Save_tbl = function(tbl, S){
  tbl.name = deparse(substitute(tbl))
  name = paste0(S,".csv")
  write.csv(tbl, paste(Exp.dr,"Qual_Data",name, sep = "/"),F)
  #append to pipeline
  cat("\nSaved",tbl.name, "as", name, file = pipetxt,append = T)
}

#_____
#LOESS_norm: Applys LOESS normalization 
#Input:   tbl - Data frame on which normalization will be performed 
#Output:  normalized dataframe
LOESS_norm = function(tbl){
  tbl.name = deparse(substitute(tbl))
  result= tbl[,!names(tbl) %in% Assays]
  for (a in Assays){
    Assay.data = na.omit(subset(tbl, select= c("Treatment","Stim","Plate","row","col",a)))
    for (p in unique(tbl$Plate)){
      #Get relevant Control data for assay and plate
        plate.data = Controls(Assay.data,a,T)
        plate.data = plate.data[plate.data$Plate == p,]
  
      #Fit the loess model on the data.
      span.val = 1
      plate.loess_model = loess(formula = get(a) ~ row + col , plate.data , span= span.val, model = T, control=loess.control(surface="direct") )
    
        newdata=data.frame(row = rep(c(1:16),each = 24),col = rep(c(1:24),16)) 
        
        newdata$extrap = predict(plate.loess_model, newdata)
        
        loess.tbl = left_join(newdata,Assay.data[Assay.data$Plate == p,c("row","col",a)],by = c("row", "col")) %>% na.omit()
        #Make data normalization step for the current plate
        med = median(loess.tbl$extrap)
        loess.tbl$result = loess.tbl[[a]] - (loess.tbl$extrap - med)
        
        loess.tbl$Plate = p
        loess.tbl %<>% select("row","col","Plate",a,"result")
        
        Assay.data = left_join(Assay.data,loess.tbl)
        Assay.data[[a]][Assay.data$Plate == p] = Assay.data$result[Assay.data$Plate == p]
        Assay.data$result = NULL
      
    }
    result = left_join(result,Assay.data,by = c("Treatment", "row", "col", "Stim", "Plate"))
  }
  #append to pipeline.txt 
  cat("\nLoess normalized:\t",tbl.name,"\n Span of",span.val, file = pipetxt,append = T)
  return(result)
}
