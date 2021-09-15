#Manuscript Plots 

#Histograms
Distributions = function(run){
    colours = c("red","blue","green")
    names(colours) = Assays
    
    #RAW data
    Raw = read.csv(paste0(Exp.dr,"/Qual_Data/RAW_data.csv")) 
    tbl = "Raw"
    plot.tbl = get(tbl)
    plots = vector("list",length(Assays))
    names(plots) = Assays
    for (a in Assays){
      (if (a == "CTG"){
        plots[[a]] = ggplot(data=plot.tbl, aes_string(a))+ 
          geom_histogram(binwidth = 100000)+ 
          geom_density(aes(y=100000* ..count..), alpha = 0.3, color = colours[a], fill = colours[a]) +
          theme_bw()+
          ggtitle(paste0(tbl," Data Distribution"))
      }
      else{
        plots[[a]] = ggplot(data=plot.tbl, aes_string(a))+ 
          geom_histogram(binwidth = 1000)+ 
          geom_density(aes(y=1000* ..count..), alpha = 0.3, color = colours[a], fill = colours[a]) +
          theme_bw()+
          ggtitle(paste0(tbl," Data Distribution"))
      })
    }
    pdf(file=paste0(Exp.dr,"/Plots/Distribution/",tbl,"_",run,".pdf"),width= 4*length(Assays),height= 2)
    grid.arrange(grobs=plots,ncol= length(Assays))
    dev.off()
    
    #Log10 data
    Log10 = read.csv(paste0(Exp.dr,"/Qual_Data/Log_Norm_data.csv"))
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
    pdf(file=paste0(Exp.dr,"/Plots/Distribution/",tbl,"_",run,".pdf"),width= 4*length(Assays),height= 2)
    grid.arrange(grobs=plots,ncol= length(Assays))
    dev.off()
    
    if (run == "LOESS"){
    #Normalized data
    Norm = read.csv(paste0(Exp.dr,"/Qual_Data/",run,"_Norm_data.csv")) 
    tbl = "Norm"
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
    pdf(file=paste0(Exp.dr,"/Plots/Distribution/",tbl,"_",run,".pdf"),width= 4*length(Assays),height= 2)
    grid.arrange(grobs=plots,ncol= length(Assays))
    dev.off()
    }
    
    #FC data
    FC = read.csv(paste0(Exp.dr,"/Qual_Data/Plate_Norm_data.csv")) 
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
    pdf(file=paste0(Exp.dr,"/Plots/Distribution/",tbl,"_",run,".pdf"),width= 4*length(Assays),height= 2)
    grid.arrange(grobs=plots,ncol= length(Assays))
    dev.off()
  
}


SSMDcrit = function (count, alpha) {
  qt((1-alpha),
     count - 1,
     ncp = sqrt(count)*0.25)*(gamma((count-1)/2)/gamma((count-2)/2))*sqrt(2/(count*(count-1)))
}



#FPL with replicates
FPL.Range = function(){
#find beta with quantile function
FPLS = data_frame(x = c(3:10,3:10), SSMD = c(SSMDcrit(3:10, 0.05),SSMDcrit(3:10,0.01)))
Statcrit0.05 = function(x){SSMDcrit(x,0.05)}
Statcrit0.01 = function(x){SSMDcrit(x,0.01)}

pl = ggplot(data.frame(x=c(3, 10)), aes(x=x)) + 
  stat_function(fun=Statcrit0.05, geom="line", aes(color = "FPL: beta2 = 0.25 \n alpha = 0.05")) + 
  stat_function(fun=Statcrit0.01, geom="line", aes(color = "FPL: beta2 = 0.25 \n alpha =  0.01")) + 
  geom_point(data = FPLS, aes(x=x,y=SSMD)) +  
  geom_label_repel(data = FPLS, aes(label = paste(round(FPLS$x,3),",",round(FPLS$SSMD,3)),x=x,y=SSMD) ,size = 2)+
  xlab("Replicate") + ylab("SSMD cut off") + 
  scale_colour_manual ("",values = c("red","blue"))+
  theme_bw()+
  theme(legend.position = c(0.85,0.95))
ggsave(paste0(Exp.dr,"/Plots/","FPL_Power",".pdf"),plot = pl, width = 5, height = 5)
}

#Biological Function Assay verification
Bio.fn = function(run){
  FC.tbl = read.csv(paste0(Exp.dr,"/Qual_Data/Plate_Norm_data.csv")) %>% filter(Treatment %in% c("NC","POS"))
  
  #CTG
  CTG.tbl = FC.tbl %>% mutate_if(is.factor, as.character)
  CTG.tbl$Treatment[CTG.tbl$Treatment == "NC"] = "No Cell"
  CTG.tbl$Treatment[CTG.tbl$Treatment == "POS"] = "Positive Control"
  CTG.tbl %<>% select(Treatment,CTG)
  write.csv(CTG.tbl, paste(Exp.dr,"Qual_Data","bio_fun_CTG.csv", sep = "/"),F)

  #LYZ.NS
  LYZ.NS.tbl = FC.tbl %>% mutate_if(is.factor, as.character)
  LYZ.NS.tbl$Stim[LYZ.NS.tbl$Treatment == "NC"] = "No Cell"
  LYZ.NS.tbl$Stim = gsub(paste(grep("^NS",unique(LYZ.NS.tbl$Stim),value = T),collapse = "|"),"Non-Stimuated",LYZ.NS.tbl$Stim)
  LYZ.NS.tbl$Stim = gsub(paste(grep("^S",unique(LYZ.NS.tbl$Stim),value = T),collapse = "|"),"CCh-Stimulated",LYZ.NS.tbl$Stim)
  LYZ.NS.tbl %<>% select(Stim,LYZ.NS)
  write.csv(LYZ.NS.tbl, paste(Exp.dr,"Qual_Data","bio_fun_LYZ.NS.csv", sep = "/"),F)
  
  #LYZ.S
  LYZ.S.tbl = FC.tbl %>% mutate_if(is.factor, as.character) %>% filter(Treatment == "NC" | Stim %in% c("NS.NS","NS.S"))
  LYZ.S.tbl$Stim[LYZ.S.tbl$Treatment == "NC"] = "No Cell"
  LYZ.S.tbl$Stim = gsub(paste(grep("NS.NS",unique(LYZ.S.tbl$Stim),value = T),collapse = "|"),"Non-Stimuated",LYZ.S.tbl$Stim)
  LYZ.S.tbl$Stim = gsub(paste(grep("NS.S",unique(LYZ.S.tbl$Stim),value = T),collapse = "|"),"CCh-Stimulated",LYZ.S.tbl$Stim)
  LYZ.S.tbl %<>% select(Stim,LYZ.S)
  write.csv(LYZ.S.tbl, paste(Exp.dr,"Qual_Data","bio_fun_LYZ.S.csv", sep = "/"),F)
  
#Plot for all plate
  plots = vector("list",3)
  
  plots[[1]]= ggplot() + 
    geom_boxplot(data = CTG.tbl,aes(x = Treatment, y = CTG , fill = Treatment),na.rm = TRUE) + 
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
  pdf(file=paste0(Exp.dr,"/Plots/","Bio_function_",run,".pdf"),width= 4 * length(Assays),height=2)
  grid.arrange(grobs=plots,ncol= 3)
  dev.off()
  
}

#Biological Function Assay verification
Bio.plate.fcn = function(run){
  FC.tbl = read.csv(paste0(Exp.dr,"/Qual_Data/Plate_Norm_data.csv")) %>% filter(Treatment %in% c("NC","POS"))
  Plates = unique(FC.tbl$Plate)
  
  #Plot for all plate
  plots = vector("list",3*length(Plates))
  
  for(p in 1: length(Plates)){
  #CTG
  CTG.tbl = FC.tbl %>% filter(Plate == Plates[p]) %>% mutate_if(is.factor, as.character)
  CTG.tbl$Treatment[CTG.tbl$Treatment == "NC"] = "No Cell"
  CTG.tbl$Treatment[CTG.tbl$Treatment == "POS"] = "Positive Control"
  CTG.tbl %<>% select(Treatment,CTG)
  write.csv(CTG.tbl, paste0(Exp.dr,"/Qual_Data/","bio_fun_CTG_",p,".csv"),F)
  
  #LYZ.NS
  LYZ.NS.tbl = FC.tbl %>% filter(Plate == Plates[p]) %>% mutate_if(is.factor, as.character)
  LYZ.NS.tbl$Stim[LYZ.NS.tbl$Treatment == "NC"] = "No Cell"
  LYZ.NS.tbl$Stim = gsub(paste(grep("^NS",unique(LYZ.NS.tbl$Stim),value = T),collapse = "|"),"Non-Stimuated",LYZ.NS.tbl$Stim)
  LYZ.NS.tbl$Stim = gsub(paste(grep("^S",unique(LYZ.NS.tbl$Stim),value = T),collapse = "|"),"CCh-Stimulated",LYZ.NS.tbl$Stim)
  LYZ.NS.tbl %<>% select(Stim,LYZ.NS)
  write.csv(LYZ.NS.tbl, paste0(Exp.dr,"/Qual_Data/","bio_fun_LYZ.NS_",p,".csv"),F)
  
  #LYZ.S
  LYZ.S.tbl = FC.tbl %>% filter(Plate == Plates[p]) %>% mutate_if(is.factor, as.character) %>% filter(Treatment == "NC" | Stim %in% c("NS.NS","NS.S"))
  LYZ.S.tbl$Stim[LYZ.S.tbl$Treatment == "NC"] = "No Cell"
  LYZ.S.tbl$Stim = gsub(paste(grep("NS.NS",unique(LYZ.S.tbl$Stim),value = T),collapse = "|"),"Non-Stimuated",LYZ.S.tbl$Stim)
  LYZ.S.tbl$Stim = gsub(paste(grep("NS.S",unique(LYZ.S.tbl$Stim),value = T),collapse = "|"),"CCh-Stimulated",LYZ.S.tbl$Stim)
  LYZ.S.tbl %<>% select(Stim,LYZ.S)
  write.csv(LYZ.S.tbl, paste0(Exp.dr,"/Qual_Data/","bio_fun_LYZ.S_",p,".csv"),F)
  
  plots[[3*(p-1)+1]]= ggplot() + 
    geom_boxplot(data = CTG.tbl,aes(x = Treatment, y = CTG , fill = Treatment),na.rm = TRUE) + 
    ggtitle(paste("Plate",p,"Viability Assay")) + labs(y = "CTG", x = "Controls") +
    theme_bw()+
    theme(legend.position = "None")+
    scale_x_discrete(limits =c("No Cell","Positive Control")) +
    scale_fill_manual(values = c("red", "blue"))
  
  plots[[3*(p-1)+2]]= ggplot() + 
    geom_boxplot(data = LYZ.NS.tbl ,aes(x = Stim, y = LYZ.NS, fill = Stim), na.rm = TRUE) + 
    ggtitle(paste("Plate",p,"Non-Stimulated Lysozyme Assay")) + labs(y = "LYZ.NS",x = "Controls") +
    theme_bw()+
    theme(legend.position = "None")+
    scale_x_discrete(limits=c("No Cell","Non-Stimuated","CCh-Stimulated")) +
    scale_fill_manual(values = c("purple", "red", "yellow"))
  
  plots[[3*(p-1)+3]]= ggplot() + 
    geom_boxplot(data = LYZ.S.tbl ,aes(x = Stim, y = LYZ.S, fill = Stim), na.rm = TRUE) + 
    ggtitle(paste("Plate",p,"CCh-Stimulated Lysozyme Assay") )+ labs(y = "LYZ.S", x = "Controls") +
    theme_bw()+
    theme(legend.position = "None")+
    scale_x_discrete(limits =c("No Cell","Non-Stimuated","CCh-Stimulated")) +
    scale_fill_manual(values = c("purple", "red", "yellow"))
  }
  
  #arrange plots and save
  pdf(file=paste0(Exp.dr,"/Plots/","Plate_Bio_function_",run,".pdf"),width= 4 * length(Assays),height=2*length(Plates))
  grid.arrange(grobs=plots,ncol= 3)
  dev.off()
  
}

#Replicate SSMD Hits
RepSSMD_plot = function(){
  #Load all SSMD Data
  SSMD.tbl = read.csv(paste0(Exp.dr,"/Analysis_Data/",Exp.folder,"_Analyzed.FC.tbl.csv")) %>% subset(!Treatment %in% c("NC","NONE","POS") & Time == "FULL",select = -X) 
  SSMD.tbl%<>% select(c("Treatment","Dose",grep("repSSMD",names(SSMD.tbl),value = T))) %>% mutate(Treatment = paste(Treatment, Dose)) %>% select(-Dose)
  names(SSMD.tbl)[names(SSMD.tbl) %in% grep("repSSMD",names(SSMD.tbl),value = T)] = sub("repSSMD\\.*","",grep("repSSMD",names(SSMD.tbl),value = T))
  
  #Load data for Rep Hits
  Hits.tbl =  read.csv(paste0(Exp.dr,"/Analysis_Data/",Exp.folder,"_Rep_hits.csv")) %>% subset(!Treatment %in% c("NC","NONE","POS")& Time == "FULL",select = -X) 
  Hits.tbl%<>% select(c("Treatment","Dose",grep("repSSMD",names(Hits.tbl),value = T))) %>% mutate(Treatment = paste(Treatment, Dose)) %>% select(-Dose)
  names(Hits.tbl)[names(Hits.tbl) %in% grep("repSSMD",names(Hits.tbl),value = T)] = sub("repSSMD\\.*","",grep("repSSMD",names(Hits.tbl),value = T))
 
  #Set limits for graph
  limits = c(min(c(SSMD.tbl$CTG,SSMD.tbl$LYZ.NS ,SSMD.tbl$LYZ.S)), max(c(SSMD.tbl$CTG,SSMD.tbl$LYZ.NS ,SSMD.tbl$LYZ.S)))*1.1
  
  #Find Treatments that are not hits
  (if (nrow(Hits.tbl) == 0){
    #append to pipeline.txt 
    cat("\nNo Replicate Hits Selected", file = pipetxt,append = T)
  }
  else{
  SSMD.tbl = setdiff (SSMD.tbl, Hits.tbl)
  })
  #Initalize list for each assay plot
  plot.col = names(SSMD.tbl)[names(SSMD.tbl) != "Treatment"]
  plots = vector("list",length(plot.col))
  
  Crit = c(SSMDcrit(8,0.05))
  colour = c("red","blue","green")
  names(colour) = plot.col
  
  for(col in 1: length(plot.col)){
    (if (col == 1){
    plots[[col]] = ggplot(data = SSMD.tbl,aes_string(x = "Treatment", y = plot.col[[col]])) + 
      geom_point(data = SSMD.tbl, color = "black") +
      geom_point(data = Hits.tbl, color = colour[col])+
      geom_hline(yintercept= 0)+
      theme_bw()+
      theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
      labs(y = "SSMD") + 
      ylim(limits) +
      ggtitle(plot.col[col])+
      geom_hline(yintercept = Crit ,linetype="dashed",color = colour[col])
    }
    else{
      plots[[col]] =  ggplot(data = SSMD.tbl, aes_string(x = "Treatment", y = plot.col[[col]])) + 
        geom_point(data = SSMD.tbl, color = "black") +
        geom_point(data = Hits.tbl, color = colour[col], size = 3)+    
        geom_hline(yintercept= 0)+
        theme_bw()+
        theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())+ 
        ylim(limits) +
        ggtitle(plot.col[col])+
        geom_hline(yintercept = Crit,linetype="dashed",color = colour[col], size = 1.5)
    })
  }
  #arrange plots and save
  pdf(file=paste0(Exp.dr,"/Plots/","RepSSMD",".pdf"),width= 4* length(plot.col),height=2)
  grid.arrange(grobs=plots,ncol= length(plot.col))
  dev.off()
}


#Z score Hits
Zscore_plot = function(){
  #Load RepSSMD hits
  Z.tbl = read.csv(paste0(Exp.dr,"/Analysis_Data/",Exp.folder,"_Rep_hits.csv")) %>% subset(!Treatment %in% c("NC","NONE","POS")& Time == "FULL",select = -X) 
  Z.tbl%<>% select(c("Treatment","Dose",grep("^Z.",names(Z.tbl),value = T))) %>% mutate(Treatment = paste(Treatment, Dose)) %>% select(-Dose)
  names(Z.tbl)[names(Z.tbl) %in% grep("^Z.",names(Z.tbl),value = T)] = sub("Z.\\.*","",grep("^Z.",names(Z.tbl),value = T))
  
  #Load Zscore hits
  Hits.tbl =  read.csv(paste0(Exp.dr,"/Analysis_Data/",Exp.folder,"_Z_Rep_hits.csv")) %>% subset(!Treatment %in% c("NC","NONE","POS")& Time == "FULL",select = -X) 
  Hits.tbl%<>% select(c("Treatment","Dose",grep("^Z.",names(Hits.tbl),value = T))) %>% mutate(Treatment = paste(Treatment, Dose)) %>% select(-Dose)
  names(Hits.tbl)[names(Hits.tbl) %in% grep("^Z.",names(Hits.tbl),value = T)] = sub("Z.\\.*","",grep("^Z.",names(Hits.tbl),value = T))
  
  #Find Treatments that are not hits
  (if (nrow(Hits.tbl) == 0){
    #append to pipeline.txt 
    cat("\nNo Z score Hits Selected", file = pipetxt,append = T)
  }
    else{
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
      names(colour) = plot.col
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
    })

}


#Dose Response Plots
Full_Dose = function(){
  #Load FC and Dose for all points
  PrimFC.tbl = read.csv(paste0(Exp.dr,"/Qual_Data/Plate_Norm_data.csv"))
  Doses = read.csv(paste0(Exp.dr,"/Analysis_Data/",Exp.folder,"_Treatment.Doses.csv")) %>% select(-X) %>% rename("1" = Dose, "2" = x2,"3" = x4,"4"= x8) %>% melt(id.vars = "Treatment") %>% rename("Dose" = value)
  #Find FC for hits
  PrimFC.tbl %<>% select("Treatment","Time","Dose","Stim","Plate", Assays) %>% full_join(Doses)
  
  #Initalize list for each assay plot
  Plots = vector("list",length(Assays))
  names(Plots) = Assays
  
  #Find Treatments
  Treatments = PrimFC.tbl$Treatment[!PrimFC.tbl$Treatment %in% c("POS","NC")] %>% unique() %>% as.character() 

  for (t in Treatments){ #for each Treatment
    
    #Find mean and variables for SE
    Dose.data = PrimFC.tbl %>% subset(Treatment %in% c(t)) 
    
    for (a in Assays){ #For each Assay
      POS.data.Early = Controls(PrimFC.tbl, a, POSCont = T) %>% select("Treatment","Dose",a) %>% mutate(Time = "EARLY")
      POS.data.Full = Controls(PrimFC.tbl, a, POSCont = T) %>% select("Treatment","Dose",a) %>% mutate(Time = "FULL")
      
      plot.data = Dose.data %>% select("Treatment","Time","Dose",a)%>% rbind(POS.data.Early,POS.data.Full)  %>%  group_by(Treatment,Time,Dose) %>% 
        summarize_at(a,funs(mean(.,na.rm = T), sd(.,na.rm =T),n()))
      Dose.data$Dose = as.numeric(sub("D","",Dose.data$Dose))
      #Select relevant data and calculate SE
      plot.data$se = plot.data$sd/sqrt(plot.data$n)
      
      plot.data = plot.data[order(plot.data$Time, plot.data$Dose),]
      plot.data$axis = 0:(length(unique(plot.data$Dose))-1)
      
      #Plot dose response curve
      Plots[[a]] =
        ggplot(data = plot.data,aes(x = axis, y = mean, colour = Time)) +
        geom_line(data=na.omit(plot.data))+
        geom_point(shape = 19, size = 1)+
        geom_errorbar(aes(ymin= mean -se, ymax= mean +se), width=.1) +
        labs(x = "Dose log2",y = "Fold Change" )+
        theme_bw()+
        theme(aspect.ratio = 3/6,text = element_text(size=15))+
        scale_x_continuous(breaks = plot.data$axis, labels = plot.data$Dose)+
        ggtitle(paste("Primary", a,":", t))
    }
    pdf(file=paste0(Exp.dr,"/Plots/Doses/",t,".pdf"),width= 6 * length(Assays),height=3)
    grid.arrange(grobs=Plots,ncol= length(Assays))
    dev.off()
  }
}


#Early Full Dose Response

Time_Dose = function(){
  #Load FC for all points
  PrimFC = read.csv(paste0(Exp.dr,"/Qual_Data/Plate_Norm_data.csv")) %>% mutate_if(is.factor,as.character)
  
  Treatments = PrimFC %>% filter(!Treatment %in% c(Cont,"NC"))
  Treatments = unique(Treatments$Treatment)
  plots = vector("list", length(Assays)*length(Treatments))
  for (t in 1:length(Treatments)){ #for each Target
    #find treatments for specific target
    Treat.data = PrimFC %>% filter(Treatment %in% c(Treatments[t],Cont))
    #Initalize list for dose response plots
    for (a in 1:length(Assays)){ # for each assay
      Col = grep(Assays[a], names(Treat.data),value = T)
      POS.data = Controls(Treat.data, Col, POSCont = T)
      Assay.data = rbind(POS.data, Treat.data[Treat.data$Treatment == Treatments[t], names(POS.data)])
      Assay.data = Assay.data %>% group_by(Dose,Treatment,Time)%>% summarize_at(Col,funs(mean(.,na.rm = T),sd(.,na.rm = T),n()))
      Assay.data %<>% mutate(se = sd/sqrt(n))
      #Plot
      plots[[((t-1)*3 + a)]] =
        ggplot() + 
        geom_point(data = Assay.data,aes(x = Time, y = mean, color = factor(Dose)),position=position_dodge(width=0.5))+
        geom_errorbar(data = Assay.data, aes(x = Time,ymin= mean-se, ymax= mean+se, color = factor(Dose)), width=.1, position=position_dodge(width=0.5))+
        stat_summary(data = Assay.data, aes(x = Time, y = mean),fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                     geom = "crossbar", width = 0.5)+
        scale_x_discrete(limits=c("NONE","EARLY","FULL")) +
        scale_colour_discrete(name = "Dose (uM)")+
        labs(y = "Fold Change\nrel. Control")+
        ggtitle(paste(Treatments[t], Assays[a]))
    }
  }
  #arrange plots and save
  pdf(file=paste(Exp.dr,"Plots","Doses","All.pdf", sep = "/"),width= 4*length(Assays),height= 2*length(plots)/3)
  grid.arrange(grobs=plots,ncol= length(Assays))
  dev.off()
}

#outlier test
POS_dist = function(run){
  #set colours for assays
  colours = c("red","blue","green")
  names(colours) = Assays
  
  #FC data
  plot.tbl = read.csv(paste0(Exp.dr,"/Qual_Data/Plate_Norm_data.csv")) %>% filter(Treatment == "POS")
  plots = vector("list",length(Assays))
  names(plots) = Assays
  for (a in Assays){
    stat =  plot.tbl %>% summarize_at(a, funs(median(.,na.rm = T),mean(.,na.rm = T)))
    plots[[a]] = ggplot(data=plot.tbl, aes_string(a))+ 
      geom_histogram(binwidth = 0.01)+ 
      geom_density(aes(y=0.01* ..count..), alpha = 0.3, color = colours[a], fill = colours[a]) +
      geom_vline(xintercept= stat$median,linetype="dashed",color = colours[a],alpha = 0.7,size = 1) +
      geom_vline(xintercept= stat$mean,color = "black",alpha = 0.7 ,size = 1) +
      theme_bw()+
      ggtitle(paste("POS data distribution","\n Mean =",round(stat$mean,3), "Median =",round(stat$median,3)) )
  }
  pdf(file=paste0(Exp.dr,"/Plots/Distribution/","POScont","_",run,".pdf"),width= 4*length(Assays),height= 2)
  grid.arrange(grobs=plots,ncol= length(Assays))
  dev.off()
  
}

#Extrapolated LOESS
LOESS_plates = function(){
  tbl = read.csv(paste0(Exp.dr,"/Qual_Data/Log_Norm_data.csv")) %>% select(-X) 
  Ident = names(tbl)[!names(tbl) %in% Assays]
  tbl =  SubLog_FC(tbl,Assays,T,tbl)
  result= tbl[,!names(tbl) %in% Assays]
  for (a in Assays){
    Assay.data = na.omit(subset(tbl, select= c("Treatment","Stim","Plate","row","col",a)))
    for (p in unique(tbl$Plate)){
        plate.data = Controls(Assay.data,a,T)
        plate.data = plate.data[plate.data$Plate == p,]
        
        pl = ggplot(Assay.data[Assay.data$Plate == p,], aes(x = col,y = -row)) +
          geom_tile(aes_string(fill = a)) +
          ggtitle(paste("Raw Data for",p)) +
          scale_fill_gradient(low = 'yellow', high = "purple", space = 'Lab')+
          scale_x_continuous(breaks = 1:24)+
          scale_y_continuous(breaks = -16:-1)+
          theme_bw()+
          theme(axis.ticks=element_blank(),panel.border=element_blank(),
                panel.grid.major=element_blank())
        ggsave(paste0(Exp.dr,"/Plots/Normalization/",a,"/Raw.",p,".pdf"),plot = pl, width = 10, height = 5)
        
      #Fit the loess model on the data.
      plate.loess_model = loess(formula = get(a) ~ row + col , plate.data , span=0.6
                                , model = T, control=loess.control(surface="direct") )
      
        POS = as.data.frame(plate.loess_model$x)
        POS$fitted = plate.loess_model$fitted
        newdata=data.frame(row = rep(c(1:16),each = 24),col = rep(c(1:24),16)) 
        
        
        Pos.fit = left_join(newdata, POS) 
        Pos.fit$fitted[is.na(Pos.fit$fitted)] = NA
        
        pl = ggplot(Pos.fit, aes(x = col,y = -row)) +
          geom_tile(aes(fill = fitted)) +
          ggtitle(paste("LOESS fit for",p)) +
          scale_fill_gradient(low = 'yellow', high = "purple", space = 'Lab')+
          scale_x_continuous(breaks = 1:24)+
          scale_y_continuous(breaks = -16:-1)+
          theme_bw()+
          theme(axis.ticks=element_blank(),panel.border=element_blank(),
                panel.grid.major=element_blank())
        
        ggsave(paste0(Exp.dr,"/Plots/Normalization/",a,"/LOESSfit.",p,".pdf"),plot = pl, width = 10, height = 5)
        
        
        newdata$extrap = predict(plate.loess_model, newdata)
        
        pl = ggplot(newdata, aes(x = col,y = -row)) +
          geom_tile(aes(fill = extrap)) +
          ggtitle(paste("Extrapolated fit for",p)) +
          scale_fill_gradient(low = 'yellow', high = "purple", space = 'Lab')+
          scale_x_continuous(breaks = 1:24)+
          scale_y_continuous(breaks = -16:-1)+
          theme_bw()+
          theme(axis.ticks=element_blank(),panel.border=element_blank(),
                panel.grid.major=element_blank())
        ggsave(paste0(Exp.dr,"/Plots/Normalization/",a,"/Extrapfit.",p,".pdf"),plot = pl, width = 10, height = 5)
        
        loess.tbl = left_join(newdata,Assay.data[Assay.data$Plate == p,c("row","col",a)],by = c("row", "col")) %>% na.omit()
        #Make data normalization step for the current plate
        med = median(loess.tbl$extrap)
        loess.tbl$result = loess.tbl[[a]] - (loess.tbl$extrap - med)
        
        loess.tbl$Plate = p
        loess.tbl %<>% select("row","col","Plate",a,"result")
        
        pl = ggplot(loess.tbl, aes(x = col,y = -row)) +
          geom_tile(aes(fill = result)) +
          ggtitle(paste("Result for",p)) +
          scale_fill_gradient(low = 'yellow', high = "purple", space = 'Lab')+
          scale_x_continuous(breaks = 1:24)+
          scale_y_continuous(breaks = -16:-1)+
          theme_bw()+
          theme(axis.ticks=element_blank(),panel.border=element_blank(),
                panel.grid.major=element_blank())
        ggsave(paste0(Exp.dr,"/Plots/Normalization/",a,"/Result.",p,".pdf"),plot = pl, width = 10, height = 5)
        
        Assay.data = left_join(Assay.data,loess.tbl)
        Assay.data[[a]][Assay.data$Plate == p] = Assay.data$result[Assay.data$Plate == p]
        Assay.data$result = NULL
    }
    result = left_join(result,Assay.data,by = c("Treatment", "row", "col", "Stim", "Plate"))
  }

  result.POS = result %>% filter(Treatment == "POS")
  plots = vector("list",length(Assays))
  names(plots) = Assays
  
  #set colours for assays
  colours = c("red","blue","green")
  names(colours) = Assays
  
  for (a in Assays){
    stat = result.POS %>% summarize_at(a, funs(median(.,na.rm = T),mean(.,na.rm = T)))
    plots[[a]] = ggplot(data=result.POS, aes_string(a))+ 
      geom_histogram(binwidth = 0.01)+ 
      geom_density(aes(y=0.01* ..count..), alpha = 0.3, color = colours[a], fill = colours[a]) +
      geom_vline(xintercept= stat$median,linetype="dashed",color = colours[a],alpha = 0.6,size = 1) +
      geom_vline(xintercept= stat$mean,color = "black",alpha = 0.8 ,size = 1) +
      theme_bw()+
      ggtitle(paste("LOESS Extrap data distribution","\n Mean =",round(stat$mean,3), "Median =",round(stat$median,3)) )
  }
  pdf(file=paste0(Exp.dr,"/Plots/Normalization/","LOESS.Extrap_POS",".pdf"),width= 4*length(Assays),height= 2)
  grid.arrange(grobs=plots,ncol= length(Assays))
  dev.off()
}

dual_flashlight= function(){
  
  SSMD.tbl = read.csv(paste0(Exp.dr,"/Analysis_Data/",Exp.folder,"_Analyzed.FC.tbl.csv")) %>% subset(!Treatment %in% c("NC","NONE","POS") & Time == "FULL",select = -X) 
  
  plots = vector("list",length(Assays))
  names(plots) = Assays
  
  for (a in Assays){
    lim_X = max(abs(SSMD.tbl[[paste0("Dose.FC.",a)]]))
    lim_Y = max(abs(SSMD.tbl[[paste0("repSSMD.",a)]]))
  plots[[a]] = ggplot(SSMD.tbl,aes_string(x =paste0("Dose.FC.",a) , y = paste0("repSSMD.",a) ))+
    #X and Y axes
    geom_hline(yintercept= 0)+
    geom_vline(xintercept= 0)+
    #points
    geom_point(shape = 16, size = 1)+
    xlim(c(-lim_X,lim_X)) +
    ylim(c(-lim_Y,lim_Y)) +
    #title plot
    ggtitle(paste("Dual Flashlight plot", a))+
    labs(x ="Fold Change Relative to Control",y = "Replicate SSMD")
  }
  
  pdf(file=paste0(Exp.dr,"/Plots/","DF_plot",".pdf"),width= 4*length(Assays),height= 3)
  grid.arrange(grobs=plots,ncol= length(Assays))
  dev.off()
}

