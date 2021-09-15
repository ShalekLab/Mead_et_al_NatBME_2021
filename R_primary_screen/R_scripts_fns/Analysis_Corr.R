#Analysis_Corr.R
#BEN MEAD, DAPHNE NOV 2017
#Doug correlation
#VERSION 1 (11/29/17)

#______
#This script will find and plot the correlation between two mice
#______
for(comb in 1: length(M.comb)){
t = "Treatment"
#Select tables with SSMDs for each mouse
mA = M.comb[[comb]][1]
mB = M.comb[[comb]][2]
mA.tbl = get(paste0("Full_",mA)) %>% separate(Plate,into=c("DM","Plate"),sep = "\\.")
mB.tbl = get(paste0("Full_",mB)) %>% separate(Plate,into=c("DM","Plate"),sep = "\\.")

#Find Columns
Cols = grep("FC", names(mA.tbl),value = T)
#Remove rows (controls, NC and NONE) and remove plate column
mA.tbl = mA.tbl[!mA.tbl[[t]] %in% c(Cont,"NC","NONE"),c(Ident,Cols)]
mB.tbl = mB.tbl[!mB.tbl[[t]] %in% c(Cont,"NC","NONE"),c(Ident,Cols)]

#Rename columns for each mouse
names(mB.tbl)[!names(mB.tbl) %in% Ident] = paste0("mB-",names(mB.tbl)[!names(mB.tbl) %in% Ident])
names(mA.tbl)[!names(mA.tbl) %in% Ident] = paste0("mA-",names(mA.tbl)[!names(mA.tbl) %in% Ident])

#Combine tables
Rep.tbl = full_join(mA.tbl,mB.tbl)
Rep.tbl %<>% na.omit

#Initialize dataframes for correlations
Corr.all = data.frame(matrix(data = NA, nrow = 1 ,ncol = length(Cols)))
names(Corr.all) = Cols
Corr = data.frame(sort(unique(Rep.tbl$Plate)),matrix(data = NA, nrow = length(unique(Rep.tbl$Plate)),ncol = length(Cols)))
names(Corr) = c("Plate",Cols)
S.Corr.all = data.frame(matrix(data = NA, nrow = 1 ,ncol = length(Cols)))
names(S.Corr.all) = Cols
S.Corr = data.frame(sort(unique(Rep.tbl$Plate)),matrix(data = NA, nrow = length(unique(Rep.tbl$Plate)),ncol = length(Cols)))
names(S.Corr) = c("Plate",Cols)

#Define colours fo each plate
colours = brewer.pal(n = length(unique(Corr$Plate)), name = "Set1")
names(colours) = unique(Rep.tbl$Plate)
Plates = unique(Rep.tbl$Plate)
plots = vector("list",length(Assays)*length(Plates))
for (c in 1:length(Cols)){
  xaxis = paste0("mA-",Cols[c])
  yaxis = paste0("mB-",Cols[c])
  #Find Pearson's correlation coefficent
  Corr.all[[Cols[c]]] = round(cor(Rep.tbl[[paste0("mA-",Cols[c])]], Rep.tbl[[paste0("mB-",Cols[c])]],method = "pearson"),3)
  #Find Spearman's correlation coefficent
  S.Corr.all[[Cols[c]]] = round(cor(Rep.tbl[[paste0("mA-",Cols[c])]], Rep.tbl[[paste0("mB-",Cols[c])]],method = "spearman"),3)
  #Plot correlation graph
  #Plot limits
  lim = c(floor(min(Rep.tbl[[xaxis]],Rep.tbl[[yaxis]])), ceiling(max(Rep.tbl[[xaxis]],Rep.tbl[[yaxis]])))
  pl = Plot_Corr(Rep.tbl,Corr.all[[Cols[c]]], S.Corr.all[[Cols[c]]] ,xaxis,yaxis, paste0(mA,".",mB,"Correlation"),paste0(mA,".",mB,"_Cor_",Cols[c]),"black",lim)
  ggsave(paste0(Exp.dr,"/Plots/Correlations/",mA,".",mB,t,"Correlation",".pdf"),plot = pl, width = 5.5, height = 5.5)
  for (p in 1: length(Plates)){#For each plate
    #Find Pearson's correlation coefficent
    Corr[[Cols[c]]][Corr$Plate == Plates[p]] = round(cor(Rep.tbl[[paste0("mA-",Cols[c])]][Rep.tbl$Plate == Plates[p]], Rep.tbl[[paste0("mB-",Cols[c])]][Rep.tbl$Plate == Plates[p]],method = "pearson"),3)
    #Find Spearman's correlation coefficent
    S.Corr[[Cols[c]]][Corr$Plate == Plates[p]] = round(cor(Rep.tbl[[paste0("mA-",Cols[c])]][Rep.tbl$Plate == Plates[p]], Rep.tbl[[paste0("mB-",Cols[c])]][Rep.tbl$Plate == Plates[p]],method = "spearman"),3)
  }
}
#Append correlation coefficents to pipeline.txt
cat("\n\nPlotted Correlation Plots for:\t",t,"between mice", mA, mB, file = pipetxt,append = T)
x = kable(Corr.all, format = "rst", row.names = FALSE)
cat("\nPearson Correlation between mice:\t",x,sep = "\n",file = pipetxt,append = T)
x = kable(Corr, format = "rst", row.names = FALSE)
cat("\nPearson Correlation per plate:",x,sep = "\n",file = pipetxt,append = T)
x = kable(S.Corr.all, format = "rst", row.names = FALSE)
cat("\nSpearman Correlation between mice:\t",x,sep = "\n",file = pipetxt,append = T)
x = kable(S.Corr, format = "rst", row.names = FALSE)
cat("\nSpearman Correlation per plate:",x,sep = "\n",file = pipetxt,append = T)

#Save coeffiecents
Save_tbl(S.Corr, paste(mA,mB,"Spearman.Corr.plate", sep = "_"))
Save_tbl(S.Corr.all, paste(mA,mB,"Spearman.Corr.all", sep = "_"))
Save_tbl(Corr, paste(mA,mB,"Pearson.Corr.plate", sep = "_"))
Save_tbl(Corr.all, paste(mA,mB,"Pearson.Corr.all", sep = "_"))

rm(c,x,p,t,xaxis,yaxis,mA.tbl,mB.tbl,Cols,mA,mB,plots,pl,Plates,lim,colours,comb,S.Corr,S.Corr.all,Corr, Corr.all,Rep.tbl)
}
