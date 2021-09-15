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
POSlists = vector("list",length(Dates))
names(POSlists) = Dates
Dlists = vector("list",length(Dates))
names(Dlists) = Dates

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
  # common start of name for assay sheets (mouse_plate)
  Headers = unlist(subset(Exp.det, Dates == Date, select = Headers))
  # start of name for compound sheets
  Comp = unlist(subset(Exp.det, Dates == Date, select = Comp))
  # start of name for Dose compound sheets
  D.Comp = unlist(subset(Exp.det, Dates == Date, select = D.Comp))
  # start of name for POS compound sheets
  POS.Comp = unlist(subset(Exp.det, Dates == Date, select = POS.Comp))
  
  #Import & save quadrant vector
  #create list of sheet names
  H.Quad = paste0(H.Quad,"_quad")
  Qlists[[Date]] = vector("list",length(Plates))
  for (i in 1:length(Plates)){
    Qlists[[Date]][[i]] = unlist(read_excel(paste(Exp.dr,"RAW_Data",File[i],sep="/"),sheet=H.Quad[i],col_names=FALSE, na = "NA"))
    Qlists[[Date]][[i]] = Qlists[[Date]][[i]][c(1,3,2,4)]
  }
  names(Qlists[[Date]]) = Plates
  
  #Import & save compound arrays
  #create list of sheet names
  Comp = paste0(Comp ,"_comp")
  #Compounds must be on 96 well plates
  Clists[[Date]]=array(0,dim = c(8,12,length(Plates)),dimnames = list(NULL,NULL,Plates))
  for (i in 1:length(Plates)){
    Clists[[Date]][,,i]=as.matrix(read_excel(paste(Exp.dr, "RAW_Data",File[i],sep="/"),sheet=Comp[i],col_names=FALSE, na = "NA"))
  }
  
  #Import & save Dose arrays
  #create list of sheet names
  D.Comp = paste0(D.Comp ,"_Dcomp")
  #Compounds must be on 96 well plates
  Dlists[[Date]]=array(0,dim = c(8,12,length(Plates)),dimnames = list(NULL,NULL,Plates))
  for (i in 1:length(Plates)){
    Dlists[[Date]][,,i]=as.matrix(read_excel(paste(Exp.dr, "RAW_Data",File[i],sep="/"),sheet=D.Comp[i],col_names=FALSE, na = "NA"))
  }
  
  #Import & save Positive control arrays
  #create list of sheet names
  POS.Comp = paste0(POS.Comp ,"_POScomp")
  #Compounds must be on 96 well plates
  POSlists[[Date]]=array(0,dim = c(8,12,length(Plates)),dimnames = list(NULL,NULL,Plates))
  for (i in 1:length(Plates)){
    POSlists[[Date]][,,i]=as.matrix(read_excel(paste(Exp.dr, "RAW_Data",File[i],sep="/"),sheet=POS.Comp[i],col_names=FALSE, na = "NA"))
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
        Alists[[Date]][[a]][,,i]=as.matrix(read_excel(paste(Exp.dr, "RAW_Data",File[i],sep="/"),sheet=Headers.Assay[i],col_names=FALSE, na = "NA"))
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

rm(Date,Plates,File,Headers,Comp,H.Quad,D.Comp,POS.Comp,first_val,i,a, Headers.Assay, File.Type)


