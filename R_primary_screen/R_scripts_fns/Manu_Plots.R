#Manuscript Plots 

#Histograms
{
  #Load Tables
  Raw.tbl = read.csv(paste0(Exp.dr,"/Qual_Data/RAW_data.csv"))
  Log10.tbl = read.csv(paste0(Exp.dr,"/Qual_Data/Log_Norm_data.csv"))
  LOESS.tbl = read.csv(paste0(Exp.dr,"/Qual_Data/LOESS_Norm_data.csv"))
  FC.tbl = read.csv(paste0(Exp.dr,"/Qual_Data/Plate_Norm_data.csv"))
  
  #For all data and each biological rep
  for(run in c(".","29","30","31")){
    #Filter data based on all/run
    Raw = Raw.tbl[grep(run, Raw.tbl$Plate),]
    Log10 = Log10.tbl[grep(run, Log10.tbl$Plate),]
    LOESS = LOESS.tbl[grep(run, LOESS.tbl$Plate),]
    FC = FC.tbl[grep(run, FC.tbl$Plate),]
    
    #set colours for assays
    colours = c("red","blue","green")
    names(colours) = Assays
    
    #RAW data
    tbl = "Raw"
    plot.tbl = get(tbl)
    plots = vector("list",length(Assays))
    names(plots) = Assays
    for (a in Assays){
      (if (a == "CTG"){
        #Histogramfor CTG assay
        plots[[a]] = ggplot(data=plot.tbl, aes_string(a))+ 
          geom_histogram(binwidth = 100000)+ 
          geom_density(aes(y=100000* ..count..), alpha = 0.3, color = colours[a], fill = colours[a]) +
          theme_bw()+
          ggtitle(paste0(tbl," Data Distribution"))
      }
      else{
        #Histogram for LYZ assays
        plots[[a]] = ggplot(data=plot.tbl, aes_string(a))+ 
          geom_histogram(binwidth = 1000)+ 
          geom_density(aes(y=1000* ..count..), alpha = 0.3, color = colours[a], fill = colours[a]) +
          theme_bw()+
          ggtitle(paste0(tbl," Data Distribution"))
      })
    }
    pdf(file=paste0(Exp.dr,"/Plots/Distribution/",tbl,"_m",run,".pdf"),width= 4*length(Assays),height= 2)
    grid.arrange(grobs=plots,ncol= length(Assays))
    dev.off()
    
    #Log10 data
    tbl = "Log10"
    plot.tbl = get(tbl)
    plots = vector("list",length(Assays))
    names(plots) = Assays
    for (a in Assays){
      plots[[a]] = ggplot(data=plot.tbl, aes_string(a))+ 
        geom_histogram(binwidth = 0.05)+ 
        geom_density(aes(y=0.05* ..count..), alpha = 0.3, color = colours[a], fill = colours[a]) +
        theme_bw()+
        ggtitle(paste0(tbl," Data Distribution"))
    }
    pdf(file=paste0(Exp.dr,"/Plots/Distribution/",tbl,"_m",run,".pdf"),width= 4*length(Assays),height= 2)
    grid.arrange(grobs=plots,ncol= length(Assays))
    dev.off()
    
    #LOESS data
    tbl = "LOESS"
    plot.tbl = get(tbl)
    plots = vector("list",length(Assays))
    names(plots) = Assays
    for (a in Assays){
      median =  plot.tbl %>% summarize_at(a, funs(median(.,na.rm = T)))
      plots[[a]] = ggplot(data=plot.tbl, aes_string(a))+ 
        geom_histogram(binwidth = 0.05)+ 
        geom_density(aes(y=0.05* ..count..), alpha = 0.3, color = colours[a], fill = colours[a]) +
        theme_bw()+
        ggtitle(paste0(tbl," Data Distribution"))
    }
    pdf(file=paste0(Exp.dr,"/Plots/Distribution/",tbl,"_m",run,".pdf"),width= 4*length(Assays),height= 2)
    grid.arrange(grobs=plots,ncol= length(Assays))
    dev.off()
    
    #FC data
    tbl = "FC"
    plot.tbl = get(tbl)
    plots = vector("list",length(Assays))
    names(plots) = Assays
    for (a in Assays){
      median =  plot.tbl %>% summarize_at(a, funs(median(.,na.rm = T)))
      plots[[a]] = ggplot(data=plot.tbl, aes_string(a))+ 
        geom_histogram(binwidth = 0.05)+ 
        geom_density(aes(y=0.05* ..count..), alpha = 0.3, color = colours[a], fill = colours[a]) +
        geom_vline(xintercept= median[[1]],linetype="dashed",color = colours[a],size = 1) +
        theme_bw()+
        ggtitle("Plate Normalized Data Distribution")
    }
    pdf(file=paste0(Exp.dr,"/Plots/Distribution/",tbl,"_m",run,".pdf"),width= 4*length(Assays),height= 2)
    grid.arrange(grobs=plots,ncol= length(Assays))
    dev.off()
  }
  
  rm(LOESS, LOESS.tbl, Log10, Log10.tbl, median, plot.tbl,Raw, Raw.tbl, a, colours, tbl, run, plots, FC, FC.tbl)
  } 

#FPL FNL Plot
{
  #Up regulated with triplicates
  count = 3 
  
  #SSMD cut off for true values to define strong or weak effects
  beta2 = 0.25
  beta1 = 3
  
  #Define variables for noncentral t distribution
  df = count - 1
  b = sqrt(count)
  k = (gamma((count-1)/2)/gamma((count-2)/2))*sqrt(2/(count*(count-1)))
  
  #Define FPL and FNL functions
  FPL = function(x) {1 - pt(x/k, df , b*beta2)}
  FNL =  function(x) {pt(x/k, df , b*beta1)}
  
  #Find interesect and important points to plot
  Error_Level = data_frame(SSMD = c(uniroot(function(x) FPL(x)-FNL(x),c(0,2))$root, 
                            uniroot(function(x) FPL(x)-0.05,c(0,2))$root,
                            uniroot(function(x) FPL(x)-0.05,c(0,2))$root,
                            uniroot(function(x) FNL(x)-0.1,c(0,2))$root,
                            uniroot(function(x) FNL(x)-0.1,c(0,2))$root
  ),
  alpha =  c(FPL(uniroot(function(x) FPL(x)-FNL(x),c(0,2))$root),
         0.05,
         FNL(uniroot(function(x) FPL(x)-0.05,c(0,2))$root),
         FPL(uniroot(function(x) FNL(x)-0.1,c(0,2))$root),
         0.1
  ),
  Error = c("FPL, FNL","FPL","FNL","FPL","FNL"))
  
  #Plot FPL and FNL functions
  pl = ggplot(data.frame(x=c(0, 2)), aes(x=x)) + 
    stat_function(fun=FPL, geom="line", aes(color = "FPL: beta2 = 0.25")) + 
    stat_function(fun=FNL, geom="line",aes(color = "FNL: beta1 = 3"))+ 
    xlab("UMVUE SSMD cutoffs") + ylab("Error Rates") + 
    geom_point(data = Error_Level, aes(x=SSMD,y=alpha)) +  
    geom_label_repel(data = Error_Level, aes(label = paste(round(Error_Level$SSMD,3),",",round(Error_Level$alpha,3)),x=SSMD,y=alpha) ,size = 2)+
    scale_colour_manual ("",values = c("red", "blue"))+
    ggtitle("Up-regulation with Triplicate") +
    theme_bw()+
    theme(legend.position = c(0.85,0.95))
  
  ggsave(paste0(Exp.dr,"/Plots/","Power",".pdf"),plot = pl, width = 5, height = 5)
  
  rm(pl, count, beta1, beta2, df, b, k, FNL, FPL)
}

#Biological Function Assay verification
{
  FC.tbl = read.csv(paste0(Exp.dr,"/Qual_Data/Plate_Norm_data.csv")) %>% filter(Treatment %in% c("NC","POS"))
  
  #Find relevant CTG controls
  CTG.tbl = FC.tbl %>% mutate_if(is.factor, as.character)
  CTG.tbl$Treatment[CTG.tbl$Treatment == "NC"] = "No Cell"
  CTG.tbl$Treatment[CTG.tbl$Treatment == "POS"] = "Positive Control"
  CTG.tbl %<>% select(Treatment,CTG)
  write.csv(CTG.tbl, paste(Exp.dr,"Qual_Data","bio_fun_CTG.csv", sep = "/"),F)

  #Find relevant LYZ.NS controls
  LYZ.NS.tbl = FC.tbl %>% mutate_if(is.factor, as.character)
  LYZ.NS.tbl$Stim[LYZ.NS.tbl$Treatment == "NC"] = "No Cell"
  LYZ.NS.tbl$Stim = gsub(paste(grep("^NS",unique(LYZ.NS.tbl$Stim),value = T),collapse = "|"),"Non-Stimuated",LYZ.NS.tbl$Stim)
  LYZ.NS.tbl$Stim = gsub(paste(grep("^S",unique(LYZ.NS.tbl$Stim),value = T),collapse = "|"),"CCh-Stimulated",LYZ.NS.tbl$Stim)
  LYZ.NS.tbl %<>% select(Stim,LYZ.NS)
  write.csv(LYZ.NS.tbl, paste(Exp.dr,"Qual_Data","bio_fun_LYZ.NS.csv", sep = "/"),F)
  
  #Find relevant LYZ.S controls
  LYZ.S.tbl = FC.tbl %>% mutate_if(is.factor, as.character) %>% filter(Treatment == "NC" | Stim %in% c("NS.NS","NS.S"))
  LYZ.S.tbl$Stim[LYZ.S.tbl$Treatment == "NC"] = "No Cell"
  LYZ.S.tbl$Stim = gsub(paste(grep("NS.NS",unique(LYZ.S.tbl$Stim),value = T),collapse = "|"),"Non-Stimuated",LYZ.S.tbl$Stim)
  LYZ.S.tbl$Stim = gsub(paste(grep("NS.S",unique(LYZ.S.tbl$Stim),value = T),collapse = "|"),"CCh-Stimulated",LYZ.S.tbl$Stim)
  LYZ.S.tbl %<>% select(Stim,LYZ.S)
  write.csv(LYZ.S.tbl, paste(Exp.dr,"Qual_Data","bio_fun_LYZ.S.csv", sep = "/"),F)
  
#Plot Bio function tables
  plots = vector("list",3)
  
  plots[[1]]= ggplot() + 
    geom_boxplot(data = CTG.tbl,aes( x = Treatment, y = CTG , fill = Treatment),na.rm = TRUE) + 
    ggtitle("Viability Assay") + labs(y = "CTG", x = "Controls") +
    theme_bw()+
    theme(legend.position = "None")+
    scale_x_discrete(limits =c("No Cell","Positive Control")) +
    scale_fill_manual(values = c("red", "blue"))
 
  plots[[2]]= ggplot() + 
    geom_boxplot(data = LYZ.NS.tbl ,aes(x = Stim, y = LYZ.NS, fill = Stim), na.rm = TRUE) + 
    ggtitle("Non-Stimulated Lysozyme Assay") + labs(y = "LYZ.NS",x = "Controls") +
    theme_bw()+
    theme(legend.position = "None")+
    scale_x_discrete(limits=c("No Cell","Non-Stimuated","CCh-Stimulated")) +
    scale_fill_manual(values = c("purple", "red", "yellow"))

  plots[[3]]= ggplot() + 
    geom_boxplot(data = LYZ.S.tbl ,aes(x = Stim, y = LYZ.S, fill = Stim), na.rm = TRUE) + 
    ggtitle("CCh-Stimulated Lysozyme Assay") + labs(y = "LYZ.S", x = "Controls") +
    theme_bw()+
    theme(legend.position = "None")+
    scale_x_discrete(limits =c("No Cell","Non-Stimuated","CCh-Stimulated")) +
    scale_fill_manual(values = c("purple", "red", "yellow"))
  
  #arrange plots and save
  pdf(file=paste0(Exp.dr,"/Plots/","Bio_function",".pdf"),width= 4 * length(Assays),height=2)
  grid.arrange(grobs=plots,ncol= 3)
  dev.off()
  
  rm(FC.tbl, CTG.tbl,LYZ.NS.tbl,LYZ.S.tbl,plots)
}

#Replicate SSMD Hits
{
  #Load all SSMD Data
  SSMD.tbl = read.csv(paste0(Exp.dr,"/Analysis_Data/TSI_Analyzed.FC.tbl.csv")) %>% subset(!Treatment %in% c("NC","NONE","POS"),select = -X) 
  SSMD.tbl%<>% select(c("Treatment","Dose",grep("_Rep",names(SSMD.tbl),value = T))) %>% mutate(Treatment = paste(Treatment, Dose)) %>% select(-Dose)
  names(SSMD.tbl)[names(SSMD.tbl) %in% grep("_Rep",names(SSMD.tbl),value = T)] = sub("SSMD.\\.*","",grep("_Rep",names(SSMD.tbl),value = T))
  names(SSMD.tbl)[names(SSMD.tbl) %in% grep("_Rep",names(SSMD.tbl),value = T)] = sub("*_Rep","",grep("_Rep",names(SSMD.tbl),value = T))
  
  #Load data for Rep Hits
  Hits.tbl =  read.csv(paste0(Exp.dr,"/Analysis_Data/TSI_Rep_hits.csv")) %>% subset(!Treatment %in% c("NC","NONE","POS"),select = -X) 
  Hits.tbl%<>% select(c("Treatment","Dose",grep("_Rep",names(Hits.tbl),value = T))) %>% mutate(Treatment = paste(Treatment, Dose)) %>% select(-Dose)
  names(Hits.tbl)[names(Hits.tbl) %in% grep("_Rep",names(Hits.tbl),value = T)] = sub("SSMD.\\.*","",grep("_Rep",names(Hits.tbl),value = T))
  names(Hits.tbl)[names(Hits.tbl) %in% grep("_Rep",names(Hits.tbl),value = T)] = sub("*_Rep","",grep("_Rep",names(Hits.tbl),value = T))
  
  #Set limits for graph
  #limits = c(min(c(SSMD.tbl$CTG,SSMD.tbl$LYZ.NS ,SSMD.tbl$LYZ.S)), max(c(SSMD.tbl$CTG,SSMD.tbl$LYZ.NS ,SSMD.tbl$LYZ.S)))
  limits = c(-3.75, max(c(SSMD.tbl$CTG,SSMD.tbl$LYZ.NS ,SSMD.tbl$LYZ.S)))
  
  #Find Treatments that are not hits
  SSMD.tbl = setdiff (SSMD.tbl, Hits.tbl)
  
  #Initalize list for each assay plot
  plot.col = names(SSMD.tbl)[names(SSMD.tbl) != "Treatment"]
  plots = vector("list",length(plot.col))
  
  #Find RepSSMD crit value and colours
  Crit = Error_Level$SSMD[1]
  colour = c("red","blue","green")
  
  #Plot Scatter plots
  for(col in 1: length(plot.col)){
    (if (col == 1){
      #For first plot keep y axis
    plots[[col]] = ggplot(data = SSMD.tbl,aes_string(x = "Treatment", y = plot.col[[col]])) + 
      geom_point(data = SSMD.tbl, color = "black") +
      geom_point(data = Hits.tbl, color = colour[col])+
      theme_bw()+
      theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
      labs(y = "SSMD") + 
      scale_y_continuous(limits = limits)+
      ggtitle(plot.col[col])+
      geom_hline(yintercept=Crit ,linetype="dashed",color = colour[col])
    }
    else{
      plots[[col]] =  ggplot(data = SSMD.tbl,aes_string(x = "Treatment", y = plot.col[[col]])) + 
        geom_point(data = SSMD.tbl, color = "black") +
        geom_point(data = Hits.tbl, color = colour[col], size = 3)+        
        theme_bw()+
        theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())+ 
        scale_y_continuous(limits = limits)+
        ggtitle(plot.col[col])+
        geom_hline(yintercept= Crit,linetype="dashed",color = colour[col], size = 1.5)
      })
  }
  #arrange plots and save
  pdf(file=paste0(Exp.dr,"/Plots/","RepSSMD",".pdf"),width= 4* length(plot.col),height=2)
  grid.arrange(grobs=plots,ncol= length(plot.col))
  dev.off()
  
  #rm(SSMD.tbl, Hits.tbl, plot.col,plots, col, Crit, colour,limits,Error_Level)
}

#Z score Hits
{
  #Load RepSSMD hits
  Z.tbl = read.csv(paste0(Exp.dr,"/Analysis_Data/TSI_Rep_hits.csv")) %>% subset(!Treatment %in% c("NC","NONE","POS"),select = -X) 
  Z.tbl%<>% select(c("Treatment","Dose",grep("^Z.",names(Z.tbl),value = T))) %>% mutate(Treatment = paste(Treatment, Dose)) %>% select(-Dose)
  names(Z.tbl)[names(Z.tbl) %in% grep("^Z.",names(Z.tbl),value = T)] = sub("Z.\\.*","",grep("^Z.",names(Z.tbl),value = T))

  #Load Zscore hits
  Hits.tbl =  read.csv(paste0(Exp.dr,"/Analysis_Data/TSI_Z_Rep_hits.csv")) %>% subset(!Treatment %in% c("NC","NONE","POS"),select = -X) 
  Hits.tbl%<>% select(c("Treatment","Dose",grep("^Z.",names(Hits.tbl),value = T))) %>% mutate(Treatment = paste(Treatment, Dose)) %>% select(-Dose)
  names(Hits.tbl)[names(Hits.tbl) %in% grep("^Z.",names(Hits.tbl),value = T)] = sub("Z.\\.*","",grep("^Z.",names(Hits.tbl),value = T))

  #Set limits
  limits = c(min(c(Z.tbl$CTG,Z.tbl$LYZ.NS ,Z.tbl$LYZ.S)), max(c(Z.tbl$CTG,Z.tbl$LYZ.NS ,Z.tbl$LYZ.S)))

  #Find Treatments that are not hits
  Z.tbl = setdiff (Z.tbl,Hits.tbl)
  
  #Initalize list for each assay plot
  plot.col = names(Z.tbl)[names(Z.tbl) != "Treatment"]
  plots = vector("list",length(plot.col))
  
  #Find Zscore crit value and colours
  Crit = qnorm(0.9)
  colour = c("red","blue","green")
  
  #Plot Scatter plots
  for(col in 1: length(plot.col)){
    (if (col == 1){
      #For first plot keep y axis
      plots[[col]] = ggplot(data = Z.tbl,aes_string(x = "Treatment", y = plot.col[[col]])) + 
        geom_point(data = Z.tbl, color = "black") +
        geom_point(data = Hits.tbl, color = colour[col])+
        theme_bw()+
        theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
        labs(y = "Z score") + 
        scale_y_continuous(limits = limits)+
        ggtitle(plot.col[col])+
        geom_hline(yintercept=Crit ,linetype="dashed",color = colour[col])
    }
    else{
      plots[[col]] =  ggplot(data = Z.tbl,aes_string(x = "Treatment", y = plot.col[[col]])) + 
        geom_point(data = Z.tbl, color = "black") +
        geom_point(data = Hits.tbl, color = colour[col], size = 3)+        
        theme_bw()+
        theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())+ 
        scale_y_continuous(limits = limits)+
        ggtitle(plot.col[col])+
        geom_hline(yintercept= Crit,linetype="dashed",color = colour[col], size = 1.5)
    })
  }
  #arrange plots and save
  pdf(file=paste0(Exp.dr,"/Plots/","Z_Rep_hits",".pdf"),width= 4 * length(plot.col),height=2)
  grid.arrange(grobs=plots,ncol= length(plot.col))
  dev.off()
  
  rm(Z.tbl, Hits.tbl, plot.col,plots, col, Crit, colour,limits)
}


#Dose Response Plots
{
  #Load Hits
  Hit.Treat = read.csv(paste0(Exp.dr,"/Analysis_Data/TSI_Z_Rep_hits.csv"))$Treatment
  #Load FC for all points
  PrimHits.tbl = read.csv(paste0(Exp.dr,"/Qual_Data/Plate_Norm_data.csv"))
  #Find FC for hits
  PrimHits.tbl %<>% subset(Treatment %in% Hit.Treat,select = c("Treatment","Dose","Plate", Assays) ) 
  
  #Save Dose Data tables
  Dose.data = PrimHits.tbl %>% rename(Rep = Plate)
  Dose.data$Rep = sub(".* ","",Dose.data$Rep)
  Dose.data$Rep = sub("\\.P.*","",Dose.data$Rep)
  Doses = as.character(unique(Dose.data$Dose))

  for (a in Assays){
  Dose.tbl = dcast(Dose.data[,c("Treatment","Rep","Dose",a)],Treatment + Rep ~ Dose, value.var = a) 
  Save_tbl(Dose.tbl, paste0(a,"_Dose_Data"))
  }
  
  #Initalize list for each assay plot
  Plots = vector("list",length(Assays))
  names(Plots) = Assays
  
  #set colours for assays
  colours = c("red","blue","green")
  names(colours) = Assays
  
  
  for (t in unique(PrimHits.tbl$Treatment)){ #for each Treatment
    #Find mean and variables for SE
    Dose.data = PrimHits.tbl %>% subset(Treatment %in% t) %>%
      group_by(Treatment,Dose) %>% 
      summarize_at(Assays,funs(mean(.,na.rm = T), sd(.,na.rm =T),n()))
    Dose.data$Dose = as.numeric(sub("D","",Dose.data$Dose))
    
    for (a in Assays){ #For each Assay
      #Select relevant data and calculate SE
      plot.data = Dose.data %>% select("Treatment","Dose", grep(a, names(Dose.data),value = T)) 
      plot.data$se = plot.data[[paste0(a,"_sd")]]/sqrt(plot.data[[paste0(a,"_n")]])
      names(plot.data)[names(plot.data) == paste0(a,"_mean")] = "meanFC"
      
      #Plot dose response curve
      Plots[[a]] =
        ggplot(data = plot.data,aes(x = Dose, y = meanFC)) +
        geom_line(color = colours[a])+
        geom_point(shape = 19, size = 1,color = colours[a])+
        geom_errorbar(aes(ymin= meanFC -se, ymax= meanFC +se), width=.1, color = colours[a]) +
        labs(x = "Dose log2",y = "Fold Change" )+
        theme_bw()+
        theme(aspect.ratio = 3/6,text = element_text(size=15))+
        scale_x_continuous(trans="log2", breaks = plot.data$Dose)+
        ggtitle(paste("Primary", a,":", t))
    }
    pdf(file=paste0(Exp.dr,"/Plots/Doses/",t,".pdf"),width= 6 * length(Assays),height=3)
    grid.arrange(grobs=Plots,ncol= length(Assays))
    dev.off()
  }
}

#Correlation plots
{
for (corr in c("Pearson.Corr.plate","Spearman.Corr.plate")){ # Plots for Pearson and Spearman correlations
  #Initailize data frame
  Corr.data = data_frame(Plate = character(),variable=character(),value=numeric(),m1=character(),m2=character())
    for (comb in M.comb){ #for each combination
    #Find relevant data and combine
    Corr.tbl = read.csv(paste(Exp.dr,"Analysis_Data",paste("TSI",paste0(comb, collapse = "_"),paste0(corr,".csv"), sep = "_"),sep = "/")) %>% select(-X)
    Corr.tbl = melt(Corr.tbl, id.vars = "Plate")
    Corr.tbl$mice = paste(comb, collapse = " ")
    Corr.data = rbind(Corr.data, Corr.tbl)
    }
  Plates = unique(Corr.data$Plate)
  Cols = unique(Corr.data$variable)
  plots = vector("list",length(Cols))
    for (c in 1: length(Cols)){ #for each Assay
      #Select Data
      plot.data = Corr.data %>% filter(variable == Cols[c])
      plot.data %<>% mutate_if(is.character, as.factor)
      #Plot mice comb vs Plate heat map
      plots[[c]] = ggplot(plot.data, aes(x = mice,y = Plate)) +
        geom_tile(aes(fill = value), color='yellow') +
        ggtitle(paste("Correlation for",Cols[c])) +
        scale_fill_gradient(low = 'yellow', high = "purple", space = 'Lab', limits = c(0,1))+
        theme_bw()+
        theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
              axis.ticks=element_blank(),panel.border=element_blank(),
              panel.grid.major=element_blank())
    }
  pdf(file=paste0(Exp.dr,"/Plots/Correlations/",corr,".pdf"),width= 4 * length(Cols),height=4)
  grid.arrange(grobs=plots,ncol= length(Cols))
  dev.off()
  
  rm(Corr.data, Corr.tbl, plot.data, c, colours, Cols, comb, corr, Plates, plots)
  }
}
