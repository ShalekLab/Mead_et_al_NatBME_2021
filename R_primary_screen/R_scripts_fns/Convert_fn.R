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
  Plist = Plists[[date]]
  #Build NC subtracted arrays
  #read Clist
  Clist = Clists[[date]]
  #list of NC subtracted dataframes for each assay
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
    B.tbls[[a]]= Arr_to_tbl(Q.B.list[[a]],Clist,a,date, Qlist, Plist)
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
Arr_to_tbl=function(Q.B.Arr,Clist,Assay,date,Qlist,Plist) {
  #Convert compound list to formatted table
  #name columns with A to L
  dimnames(Clist)[2]=list(LETTERS[1:12])
  #Create Identifying columns
  Clist= Clist %>%
    as.data.frame.table() %>%
    mutate(Index=paste(Var1,Var2)) %>%
    mutate(ID=paste(Var1,Var2,Var3))
  #Rename column
  names(Clist)[names(Clist) == "Freq"] ="Treatment"
  #Select relevant colums
  Clist =  Clist[,c("Treatment","Index","ID")]
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
    row1 = c(Qlist[[p]][c(1,2)], Plist[[p]][c(1,2)])
    row2 = c(Qlist[[p]][c(3,4)], Plist[[p]][c(3,4)])
    col1 = c(Qlist[[p]][c(1,3)], Plist[[p]][c(1,3)])
    col2 = c(Qlist[[p]][c(2,4)], Plist[[p]][c(2,4)])
    #Find 384 row and col number
    B.tbl$row[B.tbl$Plate == p & B.tbl$Quad %in% row1] = (2*B.tbl$row[B.tbl$Plate == p & B.tbl$Quad %in% row1])-1
    B.tbl$row[B.tbl$Plate == p & B.tbl$Quad %in% row2] = 2*B.tbl$row[B.tbl$Plate == p & B.tbl$Quad %in% row2]
    B.tbl$col[B.tbl$Plate == p & B.tbl$Quad %in% col1] = (2*B.tbl$col[B.tbl$Plate == p & B.tbl$Quad %in% col1])-1
    B.tbl$col[B.tbl$Plate == p & B.tbl$Quad %in% col2] = 2*B.tbl$col[B.tbl$Plate == p & B.tbl$Quad %in% col2]
  }

  #rename & format tabl3
       #Changes made here
  B.tbl=left_join(Clist,B.tbl,by = "ID")
  B.tbl$ID = NULL
  for (p in names(Qlist)){ #for each plate
    for (q in 1:4){
      B.tbl$Quad[B.tbl$Plate == p & B.tbl$Treatment %in% unique(grep("POS",B.tbl$Treatment, value = T)) & B.tbl$Quad == Qlist[[p]][q]] = Plist[[p]][q]
    }}
  B.tbl=mutate(B.tbl,Plate=paste(date,Plate))
  B.tbl = separate(B.tbl, Quad, c("Dose","Cell.Type","Stim"), "[ ]")
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
      mtbl = left_join(mtbl,Arr.list[[a]],by = c("Treatment", "Index", "row", "col", "Dose", "Cell.Type", "Stim", "Plate"))
    }
    #remove empty Control rows
    mtbl=mtbl[c(names(mtbl)[!names(mtbl) %in% Assays],Assays)]
    return(mtbl)
}


#Save_tbl: Saves CSV of table
#Input:   tbl - table to save
#         S - name to save file as
#Output: NONE - Saves .csv file in Data folder
Save_tbl = function(tbl, S){
  tbl.name = deparse(substitute(tbl))
  name = paste0(S,".csv")
  write.csv(tbl, paste(Exp.dr,"Import_Data",name, sep = "/"),F)
  #append to pipeline
  cat("\n\tSaved",tbl.name, "as", name, "in Qual_Data", file = pipetxt,append = T)
}

