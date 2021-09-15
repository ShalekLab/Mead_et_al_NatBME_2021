#Import.R
#BEN MEAD, DAPHNE SZE JUNE 2017
#VERSION 2.0 (6/22/17)
#______
#This script is for defining file and sheet names, and importing the excel files into respective arrays
#______
# Initialize lists for Quadrant list, compound list, Assay values

Qlists = vector("list",length(Dates))
names(Qlists) = Dates
Clists = vector("list",length(Dates))
names(Clists) = Dates
Alists = vector("list",length(Dates))
names(Alists) = Dates
Plists = vector("list",length(Dates))
names(Plists) = Dates

#Experiments
#   Date = Date of experiment (used for naming) MD
for (Date in as.character(Dates)){
  # list of plate names (Mouse.plate)
  Plates = unlist(subset(Exp.det, Dates == Date, select = Plates))
  # file type (csv or xlsx)
  File.Type = unlist(subset(Exp.det, Dates == Date, select = File_Type))
  # filename (to access file)
  File = unlist(subset(Exp.det, Dates == Date, select = File))
  # common start of name for quadrant sheets 
  H.Quad = unlist(subset(Exp.det, Dates == Date, select = H.Quad))
  # common start of name for all POS control quadrant sheets 
  H.POS = unlist(subset(Exp.det, Dates == Date, select = H.POS))
  # common start of name for assay sheets (mouse_plate)
  Headers = unlist(subset(Exp.det, Dates == Date, select = Headers))
  # start of name for compund sheets
  Comp = unlist(subset(Exp.det, Dates == Date, select = Comp))
  
  #Import & save quadrant vector
  #create list of sheet names
  H.Quad = paste0(H.Quad,"_quad")
  Qlists[[Date]] = vector("list",length(Plates))
  H.POS = paste0(H.POS,"_POS")
  Plists[[Date]] = vector("list",length(Plates))
  for (i in 1:length(Plates)){
    Qlists[[Date]][[i]] = unlist(read_excel(paste(Exp.dr,"RAW_Data",File[i],sep="/"),sheet=H.Quad[i],col_names=FALSE))
    Qlists[[Date]][[i]] = Qlists[[Date]][[i]][c(1,3,2,4)]
    Plists[[Date]][[i]] = unlist(read_excel(paste(Exp.dr,"RAW_Data",File[i],sep="/"),sheet=H.POS[i],col_names=FALSE))
    Plists[[Date]][[i]] = Plists[[Date]][[i]][c(1,3,2,4)]
  }
  names(Qlists[[Date]]) = Plates
  names(Plists[[Date]]) = Plates
  
  #Import & save compound arrays
  #create list of sheet names
  Comp = paste0(Comp ,"_comp")
  #Compounds must be on 96 well plates
  Clists[[Date]]=array(0,dim = c(8,12,length(Plates)),dimnames = list(NULL,NULL,Plates))
  for (i in 1:length(Plates)){
    Clists[[Date]][,,i]=as.matrix(read_excel(paste(Exp.dr, "RAW_Data",File[i],sep="/"),sheet=Comp[i],col_names=FALSE))
  }
  
  #Import & save Assay arrays
  Alists[[Date]] = vector("list",length(Assays))
  names(Alists[[Date]]) = Assays
  for (a in Assays){
    # Assays must be on 384 well plates
    Alists[[Date]][[a]] = array(0,dim = c(16,24,length(Plates)),dimnames = dimnames(Clists[[Date]]))
    first_val = vector()
    Headers.Assay =  paste(Headers,a, sep = "_")
    for (i in 1:length(Plates)){
      if(File.Type[i] == "csv"){
      Vals = read.csv(paste(Exp.dr, "RAW_Data",paste0(a,".Vals"),paste0(Headers[i],".csv"),sep="/"),header = F)
      Vals %<>% separate(V1, into = c("row","col"), sep = 1)
      Vals = Vals[order(Vals$row,as.numeric(Vals$col)),]
      Alists[[Date]][[a]][,,i]=matrix(Vals$V2, nrow = 16, ncol = 24, byrow = T)
      first_val = c(Alists[[Date]][[a]][1,1,i],first_val)
      }
      if(File.Type[i] == "xlsx"){
        Alists[[Date]][[a]][,,i]=as.matrix(read_excel(paste(Exp.dr, "RAW_Data",File[i],sep="/"),sheet=Headers.Assay[i],col_names=FALSE))
        first_val = c(Alists[[Date]][[a]][1,1,i],first_val)
      }
    }
    if(anyDuplicated(first_val) > 0){
      dup = paste(Date,Plates[first_val %in% first_val[duplicated(first_val)]],sep ="_")
      cat("\nDuplicated plates:\t",dup, file = pipetxt,append = T, sep = "\t")
      rm(dup)
    }
  }
  #append to pipeline.txt
  (if (length(grep("Imported", readLines(pipetxt))) == pipe.ran){
    cat("\nImported:\n\t",unique(File), file = pipetxt,append = T)
  }
    else{
      cat("\n\t", unique(File), file = pipetxt,append = T)
    })
  
}

rm(Date,Plates,File,Headers,Comp,H.Quad,H.POS,first_val,i,a, Vals, Headers.Assay, Exp.det)
