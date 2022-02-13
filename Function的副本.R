
#ROC survival function
ROCsurvival =  function(training_dataset,testing_dataset,file){
  #ROC
  T1 <- training_dataset %>% dplyr::rename("censor" = "OS","time" = "OStime")
  T2 <- testing_dataset %>% dplyr::rename("censor" = "OS","time" = "OStime")
  T3 <- bind_rows(T1,T2) 
  
  #ROC survival analysis
  {
    sur365 <- survivalROC(Stime=T1$time,
                          status=T1$censor,    
                          marker = T1$Score,     
                          predict.time = 365,
                          span = 0.25*nrow(T1)^(-0.20)) 
    
    sur1085 <- survivalROC(Stime=T1$time,
                           status=T1$censor,    
                           marker = T1$Score,     
                           predict.time = 1085,
                           span = 0.25*nrow(T1)^(-0.20))
    
    sur1825 <- survivalROC(Stime=T1$time,
                           status=T1$censor,    
                           marker = T1$Score,     
                           predict.time = 1825,
                           span = 0.25*nrow(T1)^(-0.20))
    
    sur365t <- survivalROC(Stime=T2$time,
                           status=T2$censor,    
                           marker = T2$Score,     
                           predict.time = 365,
                           span = 0.25*nrow(T1)^(-0.20)) 
    
    sur1085t <- survivalROC(Stime=T2$time,
                            status=T2$censor,    
                            marker = T2$Score,     
                            predict.time = 1085,
                            span = 0.25*nrow(T1)^(-0.20))
    
    sur1825t <- survivalROC(Stime=T2$time,
                            status=T2$censor,    
                            marker = T2$Score,     
                            predict.time = 1825,
                            span = 0.25*nrow(T1)^(-0.20))
    
    sur365a <- survivalROC(Stime=T3$time,
                           status=T3$censor,    
                           marker = T3$Score,     
                           predict.time = 365,
                           span = 0.25*nrow(T1)^(-0.20)) 
    
    sur1085a <- survivalROC(Stime=T3$time,
                            status=T3$censor,    
                            marker = T3$Score,     
                            predict.time = 1085,
                            span = 0.25*nrow(T1)^(-0.20))
    
    sur1825a <- survivalROC(Stime=T3$time,
                            status=T3$censor,    
                            marker = T3$Score,     
                            predict.time = 1825,
                            span = 0.25*nrow(T1)^(-0.20))
    
    png(filename = here::here(file,"Plot","training_ROCsurvival.png"),res = 300,width = 6,height = 6,units = "in")
    plot(sur365$FP, sur365$TP,
         type="l",col=ggsci::pal_nejm()(1), 
         lwd = 3,
         xlim=c(0,1), ylim=c(0,1),   
         xlab= "1 - Specificity", ##连接
         ylab="Sensitivity",
         main="Training Data Set")## \n换行符
    lines(sur1085$FP, sur1085$TP, type="l",lwd = 3,
          col=ggsci::pal_nejm()(2)[2],xlim=c(0,1), ylim=c(0,1))
    lines(sur1825$FP, sur1825$TP, type="l",lwd = 3,
          col=ggsci::pal_nejm()(3)[3],xlim=c(0,1), ylim=c(0,1))
    abline(0,1,col="gray",lty=2,
           lwd = 2)##线条颜色
    legend(0.4,0.2,c(paste("AUC of 1 year survival =",round(sur365$AUC,3)),
                     paste("AUC of 3 year survival =",round(sur1085$AUC,3)),
                     paste("AUC of 5 year survival =",round(sur1825$AUC,3))),
           x.intersp=1, y.intersp=0.8,
           lty= 1 ,lwd= 2,col=ggsci::pal_nejm()(3),
           bty = "n",# bty框的类型
           seg.len=1,cex=1)
    dev.off()
    
    png(filename = here::here(file,"Plot","testing_ROCsurvival.png"),res = 300,width = 6,height = 6,units = "in")
    plot(sur365t$FP, sur365t$TP,
         type="l",col=ggsci::pal_nejm()(1), 
         lwd = 3,
         xlim=c(0,1), ylim=c(0,1),   
         xlab= "1 - Specificity", ##连接
         ylab="Sensitivity",
         main="Testing Data Set")## \n换行符
    lines(sur1085t$FP, sur1085t$TP, type="l",lwd = 3,
          col=ggsci::pal_nejm()(2)[2],xlim=c(0,1), ylim=c(0,1))
    lines(sur1825t$FP, sur1825t$TP, type="l",lwd = 3,
          col=ggsci::pal_nejm()(3)[3],xlim=c(0,1), ylim=c(0,1))
    abline(0,1,col="gray",lty=2,
           lwd = 2)##线条颜色
    legend(0.4,0.2,c(paste("AUC of 1 year survival =",round(sur365t$AUC,3)),
                     paste("AUC of 3 year survival =",round(sur1085t$AUC,3)),
                     paste("AUC of 5 year survival =",round(sur1825t$AUC,3))),
           x.intersp=1, y.intersp=0.8,
           lty= 1 ,lwd= 2,col=ggsci::pal_nejm()(3),
           bty = "n",# bty框的类型
           seg.len=1,cex=1)
    dev.off()
    
    png(filename = here::here(file,"Plot","All_ROCsurvival.png"),res = 300,width = 6,height = 6,units = "in")
    plot(sur365a$FP, sur365a$TP,
         type="l",col=ggsci::pal_nejm()(1), 
         lwd = 3,
         xlim=c(0,1), ylim=c(0,1),   
         xlab= "1 - Specificity", ##连接
         ylab="Sensitivity",
         main="All Data Set")## \n换行符
    lines(sur1085a$FP, sur1085a$TP, type="l",lwd = 3,
          col=ggsci::pal_nejm()(2)[2],xlim=c(0,1), ylim=c(0,1))
    lines(sur1825a$FP, sur1825a$TP, type="l",lwd = 3,
          col=ggsci::pal_nejm()(3)[3],xlim=c(0,1), ylim=c(0,1))
    abline(0,1,col="gray",lty=2,
           lwd = 2)##线条颜色
    legend(0.4,0.2,c(paste("AUC of 1 year survival =",round(sur365a$AUC,3)),
                     paste("AUC of 3 year survival =",round(sur1085a$AUC,3)),
                     paste("AUC of 5 year survival =",round(sur1825a$AUC,3))),
           x.intersp=1, y.intersp=0.8,
           lty= 1 ,lwd= 2,col=ggsci::pal_nejm()(3),
           bty = "n",# bty框的类型
           seg.len=1,cex=1)
    dev.off()
  }
}

risk_survival = function(data,title,genelist,file,name){
  #NEED TIME CENSOR and GENELIST
  list <- vector(mode = "list",length = 3)
  
  data <- data %>%
    dplyr::rename("time" = "OStime",
                  "censor" = "OS")
  
  model <- coxph(Surv(time,censor) ~ Score, data = data)
  phenotype <- data[,c("time","censor")]
  riskscore <- predict(model,data) 
  
  data$riskscore <- riskscore
  data$riskgroup <- ifelse(data$riskscore > median(data$riskscore),"high","low")
  
  risk_df <- data.frame(Patients = 1:length(riskscore),
                        RiskScore=as.numeric(sort(riskscore))) 
  risk_df$group <- ifelse(risk_df$RiskScore > median(risk_df$RiskScore),'High Risk','Low Risk')
  
  survival_df <- data.frame(Patients = 1:length(riskscore),
                            Days = phenotype[names(sort(riskscore)),"time"],
                            Status = phenotype[names(sort(riskscore)),"censor"])
  
  survival_df$Status <- ifelse(survival_df$Status==0,'Alive','Death')
  
  expr_df <- data[names(sort(riskscore)),which(colnames(data) %in% genelist)]
  
  #plot point 
  plot_point <- ggplot(risk_df,aes(x=Patients,y=RiskScore,color = group))+
    geom_point() + 
    scale_color_manual(values = ggsci::pal_jco()(2)) +
    labs(fill = "Group",y = "Risk Score",title = title) +
    theme_classic() +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) 
  
  #survival plot point 
  plot_survival <- ggplot(survival_df,aes(x=Patients,y=Days)) + 
    scale_color_manual(values = ggsci::pal_jco()(2)) +
    geom_point(aes(color = Status)) +
    theme_classic() +
    theme(legend.title = element_blank()) 
  
  #heatmap
  mycolors <- colorRampPalette(c("black", "green", "red"), bias = 1.2)(100)
  heatmap_df <- t(scale(expr_df[,genelist]))
  heatmap_df[heatmap_df > 1] = 1
  heatmap_df[heatmap_df < -1] = -1
  
  plot_heatmap <- pheatmap(heatmap_df,col= mycolors,show_colnames = F,cluster_cols = F)
  
  list[[1]] <- plot_point
  list[[2]] <- plot_survival
  list[[3]] <- plot_heatmap
  
  ggsave(plot = list[[1]],filename = here::here(file,"Plot",paste0(name,"_Risk.tiff")),
         compression = "lzw",device = "tiff",
         width = 5,height = 4,dpi = 300)
  
  ggsave(plot = list[[2]],filename = here::here(file,"Plot",paste0(name,"_Survival.tiff")),
         compression = "lzw",device = "tiff",
         width = 5,height = 4,dpi = 300)
  
  ggsave(plot = list[[3]],filename = here::here(file,"Plot",paste0(name,"_heatmap.tiff")),
         compression = "lzw",device = "tiff",
         width = 5,height = 3,dpi = 300)
}
