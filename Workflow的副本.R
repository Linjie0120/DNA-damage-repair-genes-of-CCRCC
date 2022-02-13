#Install packages
#install.packages("Package")
#BiocManager::install("Pakcages")
library(tidyverse)
library(tidymodels)
library(openxlsx)
library(survival)
library(survminer)
library(glmnet)
library(survivalROC)
library(clusterProfiler)
library(enrichplot)
library(limma)
library(rms)
library(MASS)
library(forestplot)
library(pheatmap)
library(openxlsx)

rm(list = ls())
#Load ccRCC data
survival_rs_repair <- openxlsx::read.xlsx(here::here("DataFiles","survival_list.xlsx"),sheet = 2)
load(file = here::here("RData","KIRC_Data.RData"))

mydf <- survival_rs_repair %>% dplyr::filter(pvalue < 0.05)

tumor_id <- KIRC_pheno_df %>% dplyr::filter(type == "Tumor" & OStime > 30) %>% dplyr::pull(TCGA_id)
KIRC_tumor_expr <- KIRC_expr %>% dplyr::select(tumor_id) %>% t() %>% as.data.frame() %>%
  dplyr::select(mydf$Symbol) %>%
  tibble::rownames_to_column(var = "TCGA_id") %>%
  dplyr::inner_join(KIRC_pheno_df, by = "TCGA_id") %>%
  dplyr::mutate(Stage =  fct_collapse(Stage,
                                      'I-II' = c("Stage I","Stage II"),
                                      'III-IV' = c("Stage III","Stage IV")))

#Split into training and testing cohort
set.seed(1234)
split <- initial_split(KIRC_tumor_expr,prop =7/10,strata = "Stage")
training <- training(split)
testing <- testing(split)

#save(training,testing,file = here::here("RData","Split.RData"))
#---------------#
#LASSO for
#Feature selection
{
  expr_matrix_train <- training[,mydf$Symbol] %>% as.matrix()
  
  survivaldata_train <- training[,c("OStime","OS")] %>% purrr::set_names("time","status") %>% as.matrix()
  
  lasso_fit <- glmnet(x = expr_matrix_train, y = survivaldata_train,family = "cox",alpha = 1 )
  png(filename = here::here("Plot","LASSO_FIT.png"),
      units = "in",res = 600,width=5,height = 5)
  plot(lasso_fit,label = F)
  dev.off()
  
  cvfit <- cv.glmnet(x = expr_matrix_train, y = survivaldata_train,family = "cox", type.measure = "C")
  png(filename = here::here("Plot","cvLASSO.png"),
      units = "in",res = 600,width=5,height = 5)
  plot(cvfit)
  dev.off()
  
  model_lasso_min <- glmnet(x = expr_matrix_train, y = survivaldata_train, alpha = 1,
                            type.measure = 'C',family = 'cox',lambda=cvfit$lambda.min)
  
  dd <- as.data.frame(model_lasso_min$beta[,1]) %>%
    tibble::rownames_to_column(var = "symbol") 
  colnames(dd)[2] <- "coef"
  dd <- dd %>%
    dplyr::filter(coef != 0) 
  
  dd1 <- dd %>%
    dplyr::mutate(type = ifelse(coef > 0, "pos","neg"))
  coef_plot <- ggplot(dd1,aes(x = fct_reorder(symbol,coef), y = coef,fill = type)) +
    geom_bar(stat="identity", position="identity") +
    scale_fill_manual(values=ggsci::pal_nejm()(2), guide=FALSE) +
    labs(x = "Gene Symbol", y = "coef") +
    coord_flip() +
    theme_classic()
  
  ggsave(coef_plot,filename = here::here("Plot","coef_plot.tiff"),
         compression = "lzw",dpi = 600,width = 5,height = 5)
  
  write.xlsx(dd,file = here::here("DataFiles","LASSO_coef.xlsx"),overwrite = T)
  
}

#Prognosis of LASSO model
{
  predict_train <- predict(model_lasso_min, newx = expr_matrix_train, s = cvfit$lambda.min) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "TCGA_id") %>%
    dplyr::mutate(TCGA_id = training$TCGA_id)
  
  names(predict_train) <- c("TCGA_id","Score")
  
  training_dataset <- training[,c("TCGA_id","OS","OStime","PFItime","PFI",dd$symbol)] %>%
    dplyr::inner_join(predict_train,by="TCGA_id")
  
  exprdata_test <- testing[,mydf$Symbol] %>% as.matrix()
  
  predict_test <- predict(model_lasso_min, newx = exprdata_test, s=cvfit$lambda.min) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "TCGA_id") %>%
    dplyr::mutate(TCGA_id = testing$TCGA_id)
  
  names(predict_test) <- c("TCGA_id","Score")
  testing_dataset <- testing[,c("TCGA_id","OStime","OS","PFItime","PFI",dd$symbol)] %>%
    dplyr::inner_join(predict_test,by="TCGA_id")
  
  save(cvfit,model_lasso_min,dd,training_dataset,testing_dataset,
       file = here::here("RData","LASSO_rs.RData"))
}

#Risk OS PFS analysis
{
  load(here::here("RData","LASSO_rs.RData"))
  
  survival_analysis_OS_PFS = function(data,title,filename){
    temp <- data[,c("OStime","OS","Score","PFItime","PFI")] %>%
      dplyr::mutate(mygroup = if_else(Score > median(Score),
                                      "High Risk","Low Risk"))
    fit <- survfit(Surv(OStime, OS) ~ mygroup, data = temp)
    fit2 <-  survfit(Surv(PFItime, PFI) ~ mygroup, data = temp)
    plot <- ggsurvplot(fit, data = temp,conf.int=F, pval=TRUE,
                       ylab="Overall Survival",
                       xlab="Survival times (DAY)",
                       title = NULL,
                       legend.title = title,
                       legend = "top",
                       legend.labs = c(paste0("High Risk(n = ",as.data.frame(table(temp$mygroup))[1,2],")"),
                                       paste0("Low Risk(n = ",as.data.frame(table(temp$mygroup))[2,2],")")),
                       palette = "nejm")
    plot2 <- ggsurvplot(fit2, data = temp,conf.int=F, pval=TRUE,
                        ylab="Progresion Free Survival",
                        xlab="Survival times (DAY)",
                        title = NULL,
                        legend.title = title,
                        legend = "top",
                        legend.labs = c(paste0("High Risk(n = ",as.data.frame(table(temp$mygroup))[1,2],")"),
                                        paste0("Low Risk(n = ",as.data.frame(table(temp$mygroup))[2,2],")")),
                        palette = "nejm")
    ggsave(plot = plot$plot,compression = "lzw",dpi = 600,width = 5,height = 5,
           filename = here::here("Plot",paste0(filename,"_OS.tiff")))
    ggsave(plot = plot2$plot,compression = "lzw",dpi = 600,width = 5,height = 5,
           filename = here::here("Plot",paste0(filename,"_PFS.tiff")))
  }
  df <- bind_rows(training_dataset,testing_dataset)
  
  survival_analysis_OS_PFS(data = training_dataset,title = "Training Dataset",filename = "Training")
  survival_analysis_OS_PFS(data = testing_dataset,title = "Testing Dataset",filename = "Testing")
  survival_analysis_OS_PFS(data = df,title = "Entire Dataset",filename = "Entire")
  
}

#Immune Marker analysis
{
  load(here::here("RData","LASSO_rs.RData"))
  df <- training_dataset %>%
    bind_rows(testing_dataset) %>%
    dplyr::select(TCGA_id,Score) %>%
    dplyr::mutate(riskgroup = if_else(Score > median(Score),"High","Low"))
  
  LUAD_df <- KIRC_expr %>%
    dplyr::select(all_of(df$TCGA_id)) %>%
    t() %>%
    as.data.frame()
  
  df3 <- LUAD_df[,c("CTLA4","CD274","PDCD1LG2","SNCA")] %>%
    tibble::rownames_to_column(var = "TCGA_id") %>%
    dplyr::inner_join(bind_rows(training_dataset,testing_dataset),by = "TCGA_id") %>%
    dplyr::mutate(riskgroup = if_else(Score > median(Score),"High","Low")) %>%
    dplyr::rename("PDL1" = "CD274",
                  "PDL2" = "PDCD1LG2",
                  "PD1" = "SNCA")
  
  genelist <- c("PDL1","PDL2","PD1","CTLA4")
  cor_plot2 <- quickcor(df3[,dd$symbol],df3[,genelist],cor.test = TRUE,
                        type = "full",method = "pearson",
                        colours = ggsci::pal_nejm()(2)) + 
    geom_raster() +
    geom_colour(colour = "white", size = 1.5) + 
    geom_mark(r = NA,sig.thres = 0.05, size = 6, colour = "black")
  
  ggsave(plot = cor_plot2,
         filename = here::here("Plot","Correlation_ImmuneMarker_Genes.tiff"),
         width = 5,height = 5,compression = "lzw",
         device = "tiff",dpi = 600)
  
  for (i in seq_along(genelist)) {
    cor_df <- tidy(cor.test(df3$Score,df3 %>% pull(genelist[i]))) %>%
      dplyr::mutate(p = if_else(p.value < 0.001,"P < 0.001",paste("P = ",as.character(round(p.value,3)))))
    
    p <- ggplot(data = df3,aes(x = !!sym(genelist[i]),y = Score)) +
      geom_point(color = "black") +
      geom_smooth(method = "lm",se = T, color = "black") +
      labs(y = "Risk Score") +
      annotate(label = paste0("R = ",round(cor_df$estimate,3), "\n ",cor_df$p),
               geom = "text",
               x=Inf, y = Inf, vjust=1, hjust=1,
               size = 5) +
      theme_classic()
    
    ggsave(plot = p,
           filename = here::here("Plot",paste0("Correlation_Score_",genelist[i],".tiff")),
           width = 5,height = 5,compression = "lzw",
           device = "tiff",dpi = 300)
  }
  
  #Immune cell infiltration
  #load data
  load(here::here("RData","pancancer_tcga_gsva.Rdata"))
  df <- training_dataset %>%
    bind_rows(testing_dataset) %>%
    dplyr::select(TCGA_id,Score) %>%
    dplyr::mutate(sample = substr(TCGA_id,1,15)) %>%
    dplyr::inner_join(tcga_gsva[,c(1,36:63)],by = "sample") %>%
    dplyr::select(-TCGA_id) %>%
    dplyr::mutate(group = ifelse(Score > median(Score),"High","Low")) %>%
    pivot_longer(cols=3:30,
                 names_to= "celltype",
                 values_to = "NES") 
  
  p1 <- ggplot(data = df, aes(x = celltype, y = NES))+
    #geom_boxplot(aes(fill = group),position = position_dodge(1),scale = "width")+
    stat_boxplot( aes(fill = group),position = position_dodge(1),
                  geom='errorbar', linetype=1, width=0.5)+  #whiskers
    geom_boxplot(aes(fill = group),position = position_dodge(1),scale = "width",
                 outlier.shape= NA) +    
    #stat_summary(fun.y=mean, geom="point", size=2) + 
    #stat_summary(fun.data = mean_se, geom = "errorbar") +
    scale_fill_manual(values = ggsci::pal_nejm()(2),name = "Group") +
    theme_classic()+
    xlab("" )+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
    stat_compare_means(aes(group=group), label = "p.signif")
  
  ggsave(plot = p1,
         filename = here::here("Plot","Immune_cell_infiltration_RiskScore.tiff"),
         width = 15,height = 6,
         device = "tiff",dpi = 600,compression = "lzw")
}

#Immune response
{
  KIRC_TIDE <- read_csv(here::here("DataFiles","TIDE_KIRC_Predict.csv")) %>%
    dplyr::rename("TCGA_id" = "Patient") %>% .[,c(1,3,4)] %>%
    dplyr::mutate(TCGA_id = substr(TCGA_id,1,15),
                  Responder = dplyr::recode(as.character(Responder),
                                            "TRUE" = "Response",
                                            "FALSE" = "No Response")) %>%
    dplyr::distinct(TCGA_id,.keep_all = T) 
  
  {
    immune_response <- training_dataset %>%
      bind_rows(testing_dataset) %>%
      dplyr::select(TCGA_id,Score) %>%
      dplyr::inner_join(KIRC_TIDE, by = "TCGA_id") %>%
      dplyr::mutate(riskgroup = if_else(Score > median(Score),"High","Low"))
    
    risk_response <- ggplot(immune_response,aes(x = Responder, y = Score)) +
      stat_boxplot(aes(fill = Responder),geom = "errorbar") +
      geom_boxplot(aes(fill = Responder)) +
      scale_fill_manual(values = ggsci::pal_nejm()(2),guide = FALSE) +
      labs(x = NULL,y = "Risk Score") +
      theme_classic() +
      stat_compare_means(aes(label=paste0("p = ",..p.format..)),
                         method = "wilcox.test",label.x = 1.5)
    
    ggsave(plot = risk_response,
           filename = here::here("Plot","Immune_Threpay_response.tiff"),
           width = 5,height = 5,compression = "lzw",
           device = "tiff",dpi = 300)
    
    
    
    cor_df <- tidy(cor.test(immune_response$Score,immune_response$TIDE)) %>%
      dplyr::mutate(p = if_else(p.value < 0.001,"P < 0.001",paste("P = ",as.character(round(p.value,3)))))
    ggplot(data = immune_response,aes(x = TIDE,y = Score)) +
      geom_point(color = "black") +
      geom_smooth(method = "lm",se = T, color = "black") +
      annotate(label = paste0("R = ",round(cor_df$estimate,3), "\n ",cor_df$p),
               geom = "text",
               x=Inf, y = Inf, vjust=1, hjust=1,
               size = 5) +
      theme_classic()
    
    per <- immune_response %>%
      dplyr::select(Responder,riskgroup) %>%
      group_by(riskgroup,Responder) %>%
      count() %>%
      dplyr::ungroup() %>%
      dplyr::group_by(riskgroup) %>%
      mutate(percent = n/sum(n),
             label_y=cumsum(percent)) %>%
      dplyr::ungroup()
    
    per
    
    chitest<-rstatix::chisq_test(x = immune_response$riskgroup,y = immune_response$Responder)
    plot_percent <- ggplot(per,aes(x = riskgroup, y = percent,
                                   fill = factor(Responder,levels = c("Response","No Response")))) +
      geom_bar(stat = "identity") +
      geom_text(aes(y = label_y,
                    label=paste0(round(percent*100),"%")), vjust= 2, colour="black", size=12) +
      scale_fill_manual(values = rev(ggsci::pal_nejm()(2))) +
      scale_y_continuous(labels = scales::percent) +
      ggtitle(paste0("Chi-Square Test : P = ",chitest$p)) +
      labs(x = NULL,y = "Percent",fill = "") +
      theme_classic() +
      theme(plot.title = element_text(size = 15,hjust = 0.5),legend.position="top")
    
    
    ggsave(plot = plot_percent,
           filename = here::here("Plot","Immune_Threpay_Percent.tiff"),
           width = 5,height = 5,compression = "lzw",
           device = "tiff",dpi = 300)
    }
}

#Risk survival
{
  ROCsurvival(training_dataset,testing_dataset,"ccRCC_Project")
  df <- bind_rows(training_dataset,testing_dataset)
  risk_survival(data = training_dataset,title = "Training Dataset",genelist = dd$symbol,
                file = name = "Training")
  risk_survival(data = testing_dataset,title = "Testing Dataset",genelist = dd$symbol,
                file = name = "Testing")
  risk_survival(data = df,title = "Entire Dataset",genelist = dd$symbol,
                file = name = "Entire")
  
}

#Muatation analysis
{
  load(here::here("RData","LASSO_rs.RData"))
  library(TCGAmutations)
  KIRC_mutation <- tcga_load(study = "KIRC")
  df <- data.frame(Tumor_Sample_Barcode = KIRC_mutation@clinical.data$Tumor_Sample_Barcode) %>%
    dplyr::mutate(TCGA_id = substr(Tumor_Sample_Barcode,1,15)) %>%
    dplyr::inner_join(bind_rows(training_dataset,testing_dataset),by = "TCGA_id") %>%
    dplyr::mutate(riskgroup = if_else(Score > median(Score),"High","Low"))
  
  high_risk <- df %>%
    dplyr::filter(riskgroup == "High") %>%
    dplyr::pull(Tumor_Sample_Barcode)
  
  low_risk <- df %>%
    dplyr::filter(riskgroup == "Low") %>%
    dplyr::pull(Tumor_Sample_Barcode)
  
  high_risk <- subsetMaf(KIRC_mutation,tsb = high_risk)
  low_risk <- subsetMaf(KIRC_mutation,tsb = low_risk)
  
  png(filename =  here::here("Plot","High_Risk_onco_EMT.png"),width = 8,height = 6,units = "in",res = 300) #can't use ggsave
  oncoplot(maf = high_risk, top = 20,titleText = "Altered in High Risk Group")
  dev.off()
  
  png(filename =  here::here("Plot","Low_Risk_onco_EMT.png"),width = 8,height = 6,units = "in",res = 300) #can't use ggsave
  oncoplot(maf = low_risk, top = 20,titleText = "Altered in Low Risk Group")
  dev.off()
  
  png(filename =  here::here("Plot","signatureGenes_onco.png"),width = 10,height = 5,units = "in",res = 300)
  oncoplot(maf = KIRC_mutation, genes = dd$symbol,draw_titv = TRUE,
           #colors = ggsci::pal_gsea()(8),
           fontSize = 1,
           #SampleNamefontSize = 1.5,
           titleFontSize = 1.7,
           legendFontSize = 1.6,
           annotationFontSize = 1.5)
  dev.off()
  
  OncogenicPathways(maf = high_risk)
  OncogenicPathways(maf = low_risk)
  
  high.sig = oncodrive(maf = high_risk, minMut = 5, pvalMethod = 'zscore')
  plotOncodrive(res = high.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 1)
  pt.vs.rt <- mafCompare(m1 = high_risk, m2 = low_risk, m1Name = 'High Risk', m2Name = 'Low Risk', minMut = 5)
  print(pt.vs.rt)
  pt.vs.rt$results$Hugo_Symbol[1:5]
  png(filename =  here::here("Plot","mafcompare_onco.png"),width = 6,height = 5,units = "in",res = 300)
  forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.01,titleSize = 1,
             lineWidth = 1)
  dev.off()
  #Survival
  png(filename =  here::here("Plot","SETD2_Mut.png"),width = 5,height = 5,units = "in",res = 300)
  mafSurvival(maf = KIRC_mutation, genes = 'SETD2', 
              time = 'OStime', 
              Status = 'OS', clinicalData = df[,c("Tumor_Sample_Barcode","OStime","OS")],
              isTCGA = FALSE)
  dev.off()
  
  png(filename =  here::here("Plot","PRKDC_Mut.png"),width = 5,height = 5,units = "in",res = 300)
  mafSurvival(maf = KIRC_mutation, genes = 'PRKDC', 
              time = 'OStime', 
              Status = 'OS', clinicalData = df[,c("Tumor_Sample_Barcode","OStime","OS")],
              isTCGA = FALSE)
  dev.off()
  
  coOncoplot(m1 = high_risk, m2 = low_risk, m1Name = 'High Risk', m2Name = 'Low Risk', 
             genes = pt.vs.rt$results$Hugo_Symbol[1:5], removeNonMutated = TRUE)
  
  tmb_highrisk <- maftools::tmb(high_risk)
  tmb_highrisk$group <- "High Risk"
  tmb_lowrisk <- maftools::tmb(low_risk)
  tmb_lowrisk$group <- "Low Risk"
  t <- bind_rows(tmb_lowrisk,tmb_highrisk)
  
  tmb_twogroup <- ggplot(t,aes(x = fct_relevel(group,"Low Risk"), y = total_perMB_log)) +
    stat_boxplot(aes(fill = group),geom = "errorbar") +
    geom_boxplot(aes(fill = group)) +
    scale_fill_manual(values = ggsci::pal_nejm()(2),guide = FALSE) +
    labs(x = NULL,y = "log(TMB)") +
    theme_classic() +
    stat_compare_means(aes(label=paste0("P = ",..p.format..)),
                       method = "t.test",label.x = 1.5)
  
  ggsave(tmb_twogroup,filename = here::here("Plot","tmb_twogroup.tiff"),
         compression = 'lzw',width = 5,height = 5)
}

#subgroup analysis 
{
  load(here::here("RData","LASSO_rs.RData"))
  load(file = here::here("TCGA_GTEX_CCLEdata","KIRC","KIRC_Data.RData"))
  
  subgroup_survival_analysis = function(data,title,filename){
    temp <- data[,c("OStime","OS","Score")] %>%
      dplyr::mutate(mygroup = if_else(Score > median(Score),
                                      "High Risk","Low Risk"))
    fit <- survfit(Surv(OStime, OS) ~ mygroup, data = temp)
    plot <- ggsurvplot(fit, data = temp,conf.int=F, pval=TRUE,
                       ylab="Overall Survival",
                       xlab="Survival times (DAY)",
                       title = NULL,
                       legend.title = title,
                       legend = "top",
                       legend.labs = c(paste0("High Risk(n = ",as.data.frame(table(temp$mygroup))[1,2],")"),
                                       paste0("Low Risk(n = ",as.data.frame(table(temp$mygroup))[2,2],")")),
                       palette = "nejm")
    ggsave(plot = plot$plot,compression = "lzw",dpi = 600,width = 5,height = 5,
           filename = here::here("Plot",
                                 paste0("subgroup_survival_analysis_",filename,".tiff")))
    return(plot)
  }
  
  df <- training_dataset %>%
    bind_rows(testing_dataset) %>%
    dplyr::inner_join(KIRC_pheno_df[,1:8], by = "TCGA_id") %>%
    dplyr::mutate(Stage =  fct_collapse(Stage,
                                        'I-II' = c("Stage I","Stage II"),
                                        'III-IV' = c("Stage III","Stage IV")),
                  age = if_else(Age > 65, ">65","<65"),
                  stage_T = fct_collapse(stage_T,
                                         "T1-T2" = c("T1","T1a","T1b","T2","T2a","T2b"),
                                         "T3-T4" = c("T3","T3b","T3a","T3c","T4")),
                  Grade = fct_collapse(Grade,
                                       "G1-G2" = c("G1","G2"),
                                       "G3-G4" = c("G3","G4"),
                                       "unknown" = "GX"),
                  stage_M = fct_collapse(stage_M,
                                         "unknown" = "MX"),
                  stage_N = fct_collapse(stage_N,
                                         "unknown" = "NX"),
                  Sex = str_to_lower(Sex),
                  group = if_else(Score > median(Score),"High","Low"))
  
  df1 <- df
  df1[df1 == "unknown"] <- NA
  
  tbl <- 
    tbl_summary(
      df1[,c("age","group","stage_T","Stage","stage_N","stage_M","Sex","Grade")],
      by = group, # split table by group
      missing = "no" # don't list missing data separately
    ) %>%
    add_n() %>% # add column with total number of non-missing observations
    add_p() %>% # test for a difference between groups
    modify_header(label = "**Variable**") %>% # update the column header
    bold_labels() 
  
  tbl %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = here::here("DataFiles","tab_remove_unknown.docx"))
  
  
  
  #prepare data
  subgroup_age <- df %>% dplyr::filter(age == ">65") 
  subgroup_age2 <- df %>% dplyr::filter(age == "<65")
  subgroup_stage <- df %>% dplyr::filter(Stage == "I-II")
  subgroup_stage2 <- df %>% dplyr::filter(Stage == "III-IV")
  subgroup_sex <- df %>% dplyr::filter(Sex == "male")
  subgroup_sex2 <- df %>% dplyr::filter(Sex == "female")
  subgroup_staget <- df %>% dplyr::filter(stage_T == "T1-T2")
  subgroup_staget2 <- df %>% dplyr::filter(stage_T == "T3-T4")
  subgroup_stagen <- df %>% dplyr::filter(stage_N == "N0")
  subgroup_stagen2 <- df %>% dplyr::filter(stage_N == "N1")
  subgroup_stagem <- df %>% dplyr::filter(stage_M == "M0")
  subgroup_stagem1 <- df %>% dplyr::filter(stage_M == "M1")
  subgroup_grade <- df %>% dplyr::filter(Grade == "G1-G2")
  subgroup_grade1 <- df %>% dplyr::filter(Grade == "G3-G4")
  
  #survival visulize
  plot1 <- subgroup_survival_analysis(data = subgroup_age, title = "Age > 65",filename = "agehigh")
  plot2 <- subgroup_survival_analysis(data = subgroup_age2, title = "Age =< 65",filename = "agelow")
  plot3 <- subgroup_survival_analysis(data = subgroup_stage, title = "Stege I-II",filename = "I-II")
  plot4 <- subgroup_survival_analysis(data = subgroup_stage2, title = "Stage III-IV",filename = "III-IV")
  plot5 <- subgroup_survival_analysis(data = subgroup_sex, title = "Male",filename = "male")
  plot6 <- subgroup_survival_analysis(data = subgroup_sex2, title = "Female",filename = "female")
  plot7 <- subgroup_survival_analysis(data = subgroup_staget, title = "T1-T2",filename = "T1-T2")
  plot8 <- subgroup_survival_analysis(data = subgroup_staget2, title = "T3-T4",filename = "T3-T4")
  plot9 <- subgroup_survival_analysis(data = subgroup_stagen, title = "N0",filename = "N0")
  plot10 <- subgroup_survival_analysis(data = subgroup_stagen2, title = "N1",filename = "N1")
  plot11 <- subgroup_survival_analysis(data = subgroup_stagem, title = "M0",filename = "M0")
  plot12 <- subgroup_survival_analysis(data = subgroup_stagem1, title = "M1",filename = "M1")
  plot13 <- subgroup_survival_analysis(data = subgroup_grade, title = "G1-G2",filename = "G1-G2")
  plot14 <- subgroup_survival_analysis(data = subgroup_grade1, title = "G3-G4",filename = "G3-G4")
}

#Risk score and DNA repair score estimated by GSVA
{
  library(GSVA)
  library(msigdbr)
  load(here::here("RData","LASSO_rs.RData"))
  kegg_db <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG") %>%
    dplyr::filter(str_detect(gs_name,pattern = "REPAIR|JOINING|HOMOLOGOUS|MICROHOMOLOGY")) %>%
    dplyr::select(gs_name,gene_symbol) %>%
    purrr::set_names("Type","Symbol") %>%
    dplyr::mutate(Type = str_replace(Type,"KEGG_",""))
  
  marker <- lapply(split(kegg_db,kegg_db$Type), function(x){
    dd = x$Symbol
    unique(dd)
  })
  
  expr <- as.matrix(KIRC_expr)
  gsva_data <- gsva(expr,marker, method = "ssgsea")
  gsva_df <- t(gsva_data) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "TCGA_id") %>%
    dplyr::inner_join(bind_rows(training_dataset,testing_dataset),by = "TCGA_id") %>%
    dplyr::mutate(group = if_else(Score > median(Score),"High","Low"))
  
  list <- colnames(gsva_df)[2:6]
  for (i in seq_along(list)) {
    p <- ggplot(gsva_df,aes(x = group, y = !!sym(list[i]))) +
      stat_boxplot(geom = "errorbar") +
      geom_boxplot(aes(fill = group)) +
      scale_fill_manual(values = ggsci::pal_nejm()(2),guide = "none") +
      labs(x = NULL,y = "GSVA Score",title = list[i]) +
      theme_classic() +
      stat_compare_means(aes(label=paste0("p = ",..p.format..)),
                         method = "wilcox.test",label.x = 1.5) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(p,filename = here::here("Plot",paste0("GSVA_Score_",list[i],".tiff")),
           compression = "lzw",dpi = 600,width = 5,height = 5)
  }
  
  #GSVA for hallmark geneset
  library(pheatmap)
  hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    .[,c(3,4)] %>%
    purrr::set_names("Set","Genes") %>%
    dplyr::mutate(Set = str_replace_all(Set,"HALLMARK_",""),
                  Set = str_replace_all(Set,"_"," "))
  
  hallmark_marker <- lapply(split(hallmark,hallmark$Set), function(x){
    dd = x$Genes
    unique(dd)
  })
  
  expr <- as.matrix(KIRC_expr)
  gsva_hallmark<- gsva(expr,hallmark_marker, method = "gsva")
  gsva_hallmark_df <- t(gsva_hallmark) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "TCGA_id") %>%
    dplyr::inner_join(bind_rows(training_dataset,testing_dataset),by = "TCGA_id") %>%
    dplyr::mutate(group = if_else(Score > median(Score),"High","Low"))
  
  heatmap_df <- gsva_hallmark_df %>% dplyr::select(TCGA_id,unique(hallmark$Set)) %>%
    tibble::column_to_rownames(var = "TCGA_id") %>%
    t() 
  annotation_df <- data.frame(group = gsva_hallmark_df$group)
  rownames(annotation_df) <- gsva_hallmark_df$TCGA_id
  annotation_df <- annotation_df %>%
    dplyr::arrange(group)
  pheatmap(heatmap_df[,rownames(annotation_df)],
           show_colnames = F,
           cluster_cols = F,
           annotation_col = annotation_df,
           color =colorRampPalette(c("blue", "white","red"))(100))
  
}

#DEGs and GO KEGG analysis
{
  load(here::here("RData","LASSO_rs.RData"))
  df <- training_dataset %>%
    bind_rows(testing_dataset) %>%
    dplyr::select(TCGA_id,Score) %>%
    dplyr::mutate(riskgroup = if_else(Score > median(Score),"High","Low"))
  
  expr <- KIRC_expr %>%
    dplyr::select(all_of(df$TCGA_id))
  
  #group list
  group_list <- factor(df$riskgroup,levels = c("Low","High"))
  
  #design
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  head(design)
  rownames(design)=colnames(expr)
  head(design)
  
  #contrast matrix
  contrast.matrix <- makeContrasts(High - Low,
                                   levels=design)
  fit <- lmFit(expr, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)   
  fit2 <- eBayes(fit2)
  allDiffgenes <- topTable(fit2,adjust="BH",number=Inf,coef=1)  %>%
    tibble::rownames_to_column(var = "Symbol")
  
  allDiffgenes %>% dplyr::filter(abs(logFC) > 1 & adj.P.Val < 0.05 )  %>% count()
  
  DEG_volcano <- ggplot(data=allDiffgenes, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(data = subset(allDiffgenes, abs(logFC) < 1 | adj.P.Val > 0.05),
               alpha=0.4, size=0.8,col = "grey")+
    geom_point(data=subset(allDiffgenes, logFC > 1 & adj.P.Val < 0.05),
               alpha=0.4, size=3,col=ggsci::pal_nejm()(2)[1])+
    geom_point(data=subset(allDiffgenes, logFC < -1 & adj.P.Val < 0.05),
               alpha=0.4, size=3,col=ggsci::pal_nejm()(2)[2])+
    geom_hline(yintercept= -log10(0.05), linetype="dashed", size = 1,lwd=0.6,alpha=0.8) + 
    geom_vline(xintercept = c(-1,1), linetype="dashed", size = 1,lwd=0.6,alpha=0.8) +
    labs(x="log2 (FoldChange)",y="-log10 (adj.P.Val)",title = "") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(colour = FALSE)
  
  ggsave(filename =  here::here("Plot","DEG_volcano.tiff"),
         plot = DEG_volcano,device = "tiff",dpi = 600,width = 5,height = 5)
  
  write.xlsx(allDiffgenes,file = here::here("DataFiles","Deg_between_two_riskgroup.xlsx"))
  
  degs <- allDiffgenes %>%
    dplyr::filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
    dplyr::pull(Symbol)
  
  degs <- bitr(degs,fromType="SYMBOL", toType = c("ENTREZID"), OrgDb="org.Hs.eg.db")
  
  go_result <- enrichGO(gene = degs$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
  kegg_result <- enrichKEGG(gene = degs$ENTREZID,pvalueCutoff = 0.05)
  
  go_plot <- dotplot(go_result,showCategory = 20) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=30)) +
    labs(title = "")
  KEGG_plot <- dotplot(kegg_result,showCategory = 20) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=30))
  
  ggsave(plot = go_plot,filename = here::here("Plot","GO_plot.tiff"),compression = "lzw",
         width = 8,height = 8,dpi = 600)
  ggsave(plot = KEGG_plot,filename = here::here("Plot","KEGG_plot.tiff"),compression = "lzw",
         width = 8,height =8,dpi = 600)
  
  #GSEA
  {
    library(msigdbr)
    c7.sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG") %>%
      .[,3:4] %>%
      purrr::set_names("term","gene")
    geneList=allDiffgenes$logFC
    names(geneList)=allDiffgenes$Symbol
    geneList=sort(geneList,decreasing = T)
    gsea_results <- GSEA(geneList, TERM2GENE=c7.sets, verbose=FALSE)
    gsea_rs <- as.data.frame(gsea_results)
    
    gsea_plot <- gseaplot2(gsea_results,geneSetID = gsea_rs$ID[1:10])
    ggsave(plot = gsea_plot,filename = here::here("Plot","KEGG_GSVA_plot.tiff"),compression = "lzw",
           width = 16,height =8,dpi = 600)
  }
  
  #GO chord
  {
    library(circlize)
    go_df <-setReadable(go_result,OrgDb = org.Hs.eg.db, keyType="ENTREZID") %>% as.data.frame()
    
    y<-go_df %>% dplyr::pull(geneID) %>% str_split("/")
    names(y) <- go_df$ID
    y <- lapply(y, "length<-",max(lengths(y)))
    yy <- as.data.frame(y)
    go_tidy <- yy %>% tidyr::pivot_longer(cols = 1:ncol(yy),names_to = "ID",values_to = "genes") %>%
      dplyr::filter(!is.na(genes)) %>%
      dplyr::mutate(ID = str_replace_all(ID,"\\.",":")) %>%
      dplyr::inner_join(go_df[,c(2,3,6)],by = "ID")
    mylist <- go_df$Description[1:5]
    t<-go_tidy %>% dplyr::filter(Description %in% mylist) %>% 
      dplyr::group_by(genes) %>% 
      dplyr::count(genes) %>% dplyr::arrange(desc(n)) %>% dplyr::filter(n > 2)
    go_chord <- go_tidy %>% dplyr::filter(genes %in% t$genes & Description %in% mylist) %>% .[,c(2,3)] %>%
      dplyr::group_by(genes,Description) %>% dplyr::count() %>%
      dplyr::mutate(Description = str_wrap(Description,20))
    
    png(filename =  here::here("Plot","GO_Chord.png"),
        width = 10,height = 10,units = "in",res = 600) 
    chordDiagram(go_chord,
                 transparency = 0.1,
                 link.lwd = 1,    # Line width
                 link.lty = 1,    # Line type
                 link.border = 1) # Border color)
    dev.off()
  }
  
  #KEGG chord
  {
    #KEGG
    KEGG_df <-setReadable(kegg_result,OrgDb = org.Hs.eg.db, keyType="ENTREZID") %>% as.data.frame()
    
    y<-KEGG_df %>% dplyr::pull(geneID) %>% str_split("/")
    names(y) <- KEGG_df$ID
    y <- lapply(y, "length<-",max(lengths(y)))
    yy <- as.data.frame(y)
    kegg_tidy <- yy %>% tidyr::pivot_longer(cols = 1:ncol(yy),names_to = "ID",values_to = "genes") %>%
      dplyr::filter(!is.na(genes)) %>%
      dplyr::mutate(ID = str_replace_all(ID,"\\.",":")) %>%
      dplyr::inner_join(KEGG_df[,c(1,2,5)],by = "ID")
    mylist <- KEGG_df$Description
    t<-kegg_tidy %>% dplyr::filter(Description %in% mylist) %>% 
      dplyr::group_by(genes) %>% 
      dplyr::count(genes) %>% dplyr::arrange(desc(n)) %>% dplyr::filter(n > 1)
    kegg_chord <- kegg_tidy %>% dplyr::filter(genes %in% t$genes & Description %in% mylist) %>% .[,c(2,3)] %>%
      dplyr::group_by(genes,Description) %>% dplyr::count() %>%
      dplyr::mutate(Description = str_wrap(Description,20))
    
    png(filename =  here::here("Plot","KEGG_Chord.png"),
        width = 8,height = 8,units = "in",res = 600) 
    chordDiagram(kegg_chord,
                 transparency = 0.1,
                 link.lwd = 1,    # Line width
                 link.lty = 1,    # Line type
                 link.border = 1) # Border color)
    dev.off()
  }
}

#Cox and Nomogram
{
  load(file = here::here("TCGA_GTEX_CCLEdata","KIRC","KIRC_Data.RData"))
  load(here::here("RData","LASSO_rs.RData"))
  
  df <- training_dataset %>%
    bind_rows(testing_dataset) %>%
    dplyr::inner_join(KIRC_pheno_df[,1:8], by = "TCGA_id") %>%
    dplyr::mutate(Stage =  fct_collapse(Stage,
                                        'I-II' = c("Stage I","Stage II"),
                                        'III-IV' = c("Stage III","Stage IV")),
                  age = if_else(Age > 65, ">65","<65"),
                  stage_T = fct_collapse(stage_T,
                                         "T1-T2" = c("T1","T1a","T1b","T2","T2a","T2b"),
                                         "T3-T4" = c("T3","T3b","T3a","T3c","T4")),
                  Grade = fct_collapse(Grade,
                                       "G1-G2" = c("G1","G2"),
                                       "G3-G4" = c("G3","G4")),
                  Sex = str_to_lower(Sex),
                  riskScore = if_else(Score > median(Score),"High","Low"),
                  riskScore = fct_relevel(riskScore,"Low"))
  
  #prepare function
  cox_rs = function(data,group){
    #temp <- data[,c("OS","OStime",group)]
    formula <- as.formula(paste('Surv(OStime, OS)~', group))
    summary_model <- summary(coxph(formula,data = data))
    temp <- as.data.frame(summary_model$coefficients) %>%
      tibble::rownames_to_column(var = "Group") %>%
      dplyr::rename(pvalue = "Pr(>|z|)",
                    se = "se(coef)",
                    exp = "exp(coef)") %>%
      mutate(pvalue = round(pvalue,digits = 6)) %>%
      mutate('Hazard ratio' = paste0(as.character(round(exp,digits = 2)),"(",
                                     round(exp-se,digits = 2),"-",round(exp+se,digits = 2),")")) %>%
      mutate(lower = exp - se,
             upper = exp + se)
    temp[1,1] <- group
    return(temp)
  }
  
  #prepare data
  group_age <- df %>% dplyr::filter(age %in% c("<65",">65")) %>% dplyr::mutate(Age = fct_relevel(age,"<65")) %>% 
    cox_rs(data = .,group = "age")
  group_stage <- df %>% dplyr::filter(Stage %in% c("I-II","III-IV")) %>% 
    dplyr::mutate(stage = droplevels(Stage)) %>%
    cox_rs(data = .,group = "Stage")
  group_sex <- df %>% dplyr::filter(Sex %in% c("female","male")) %>% 
    cox_rs(data = .,group = "Sex")
  group_staget <- df %>% dplyr::filter(stage_T %in% c("T1-T2","T3-T4")) %>% 
    dplyr::mutate(stage_T = droplevels(stage_T)) %>%
    cox_rs(data = .,group = "stage_T")
  group_stagen <- df %>% dplyr::filter(stage_N %in% c("N0","N1")) %>% 
    dplyr::mutate(stage_N = droplevels(as.factor(stage_N))) %>%
    cox_rs(data = .,group = "stage_N")
  group_stagem <- df %>% dplyr::filter(stage_M %in% c("M0","M1")) %>% 
    dplyr::mutate(stage_M = droplevels(as.factor(stage_M))) %>%
    cox_rs(data = .,group = "stage_M")
  group_grade <- df %>% dplyr::filter(Grade %in% c("G1-G2","G3-G4")) %>% 
    dplyr::mutate(Grade = droplevels(as.factor(Grade))) %>%
    cox_rs(data = .,group = "Grade")
  group_riskscore <- df %>% cox_rs(data = .,group = "riskScore")
  
  #merge the result
  univariate_rs <- bind_rows(group_age,group_sex,group_stage,group_staget,group_stagen,group_stagem,group_grade,group_riskscore) %>%
    dplyr::mutate(Pvalue = if_else(pvalue < 0.001,"< 0.001",as.character(round(pvalue,4)))) 
  
  univariate_rs_df <- dplyr::select(univariate_rs,Group,Pvalue,`Hazard ratio`)
  
  labeltext_univariate_rs <- as.data.frame(lapply(univariate_rs_df,as.character)) %>%
    dplyr::mutate(Group = dplyr::recode(Group,"Sex" = "Sex",
                                        "Stage" = "Stage",
                                        "age" = "Age",
                                        "stage_T" = "T Stage",
                                        "stage_N" = "N Stage",
                                        "stage_M"  ="M Stage",
                                        "riskScore" = "Risk Score")) %>%
    add_row(Group = c(" "),
            Pvalue = c("P Value"),
            'Hazard.ratio' = c("Hazard Ratio"),
            .before = 1)
  
  png(filename =  here::here("Plot","univariate_rs.png"),width = 6,height = 6,units = "in",res = 300) #can't use ggsave
  p1 <- forestplot::forestplot(labeltext = labeltext_univariate_rs,
                               mean = c(NA,univariate_rs[,"exp"]),
                               lower = c(NA,univariate_rs[,"lower"]),
                               upper = c(NA,univariate_rs[,"upper"]),
                               zero = 1,
                               boxsize = 0.2,
                               lwd.zero = 2,
                               lwd.xaxis = 2,
                               lty.ci = "solid",
                               ci.vertices=TRUE, ci.vertices.height = 0.2,
                               col = fpColors(box="black", lines="black", zero = "gray50"),
                               graph.pos = 4,
                               #txt_gp = fpTxtGp(label = gpar(fontfamily = "Times")),
                               title = "Univariate COX regression")
  dev.off()
  
  #---------#
  #Multivariate cox
  mydf <- df %>% dplyr::filter(Stage %in% c("I-II","III-IV")) %>% 
    dplyr::filter(age %in% c(">65","<65")) %>%
    dplyr::filter(stage_T %in% c("T1-T2","T3-T4")) %>% 
    dplyr::filter(Grade %in% c("G1-G2","G3-G4")) %>%
    dplyr::filter(stage_N %in% c("N0","N1")) %>% 
    dplyr::filter(stage_M %in% c("M0","M1")) %>% 
    dplyr::mutate(across(where(is.character),as.factor)) %>%
    dplyr::mutate(across(where(is.factor),droplevels)) %>%
    dplyr::mutate(#Stage = droplevels(Stage),
      #stage_N = droplevels(stage_N),
      #stage_T = droplevels(stage_T),
      #stage_M = droplevels(stage_M),
      riskScore = if_else(Score > median(Score),"High","Low"),
      riskScore = fct_relevel(riskScore,"Low"))
  
  res.cox <- coxph(Surv(OStime,OS)~age+Grade+Stage+stage_T+stage_M+stage_N+riskScore,data = mydf)
  step <- stepAIC(res.cox,direction = "both")
  best.cox <- coxph(Surv(OStime, OS) ~ age + Grade + Stage + stage_M + riskScore,data = mydf)
  
  best_rs <- cox_rs(data = mydf,group = "age + Grade + Stage + stage_M + riskScore") %>%
    dplyr::mutate(Pvalue = if_else(pvalue < 0.001,"< 0.001",as.character(round(pvalue,4)))) 
  
  multivariate_rs_df <- dplyr::select(best_rs,Group,Pvalue,`Hazard ratio`)
  
  labeltext_multivariate_rs <- as.data.frame(lapply(multivariate_rs_df,as.character)) %>%
    dplyr::mutate(Group = dplyr::recode(Group,
                                        "age + Grade + Stage + stage_M + riskScore" = "Age",
                                        "GradeG3-G4" = "Grade",
                                        "StageIII-IV" = "Stage",
                                        "stage_MM1" = "M Stage",
                                        "riskScoreHigh" = "Risk Score")) %>%
    add_row(Group = c(" "),
            Pvalue = c("P Value"),
            'Hazard.ratio' = c("Hazard Ratio"),
            .before = 1)
  
  png(filename =  here::here("Plot","multivariate_rs.png"),width = 6,height = 6,units = "in",res = 300) #can't use ggsave
  p1 <- forestplot::forestplot(labeltext = labeltext_multivariate_rs,
                               mean = c(NA,best_rs[,"exp"]),
                               lower = c(NA,best_rs[,"lower"]),
                               upper = c(NA,best_rs[,"upper"]),
                               zero = 1,
                               boxsize = 0.2,
                               lwd.zero = 2,
                               lwd.xaxis = 2,
                               lty.ci = "solid",
                               ci.vertices=TRUE, ci.vertices.height = 0.2,
                               col = fpColors(box="black", lines="black", zero = "gray50"),
                               graph.pos = 4,
                               #txt_gp = fpTxtGp(label = gpar(fontfamily = "Times")),
                               title = "Multivariate COX regression")
  dev.off()
  
  #Nomogram
  library(rms)
  
  mydf <- df %>% dplyr::filter(Stage %in% c("I-II","III-IV")) %>% 
    dplyr::filter(age %in% c(">65","<65")) %>%
    dplyr::filter(stage_T %in% c("T1-T2","T3-T4")) %>% 
    dplyr::filter(Grade %in% c("G1-G2","G3-G4")) %>%
    dplyr::filter(stage_N %in% c("N0","N1")) %>% 
    dplyr::filter(stage_M %in% c("M0","M1")) %>% 
    dplyr::mutate(across(where(is.character),as.factor)) %>%
    dplyr::mutate(across(where(is.factor),droplevels)) %>%
    dplyr::mutate(#Stage = droplevels(Stage),
      #stage_N = droplevels(stage_N),
      #stage_T = droplevels(stage_T),
      #stage_M = droplevels(stage_M),
      riskScore = if_else(Score > median(Score),"High","Low"),
      riskScore = fct_relevel(riskScore,"Low"))
  
  nomogram_data <- mydf
  
  label(nomogram_data$Score) <- "Risk Score"
  label(nomogram_data$stage_M) <- "M Stage"
  label(nomogram_data$stage) <- "Stage"
  label(nomogram_data$age) <- "Age"
  label(nomogram_data$Score) <- "Risk Score"
  label(nomogram_data$Grade) <- "Grade"
  
  ddist <- datadist(nomogram_data); options(datadist='ddist') #necessary step 
  
  a <- cph(Surv(OStime, OS) ~ age + Grade + Stage + stage_M + Score,data=nomogram_data,
           surv=TRUE,x=TRUE,y=TRUE)
  summary(a)
  surv <- Survival(a)
  nom <- nomogram(a, fun=list(function(x) surv(365, x),
                              function(x) surv(1095, x),
                              function(x) surv(1825, x)),
                  funlabel=c("Probability of \n1 year survival", 
                             "Probability of \n3 years survival",
                             "Probability of \n5 years survival"), lp=F)
  
  png(filename = here::here("Plot","Nomogram.png"),width = 15,height = 10,res = 300,units = "in")
  plot(nom, xfrac=.2, 
       total.points.label="Sum of all points", 
       cex.axis = 1.05,
       #force.label = TRUE,
       tcl = 0.8,
       lmgp = 0.1,
       vnames="labels",
       col.grid=gray(c(0.85,0.95)))
  dev.off()
  
  #----#
  source(file = here::here("Code","function_cal.R"))
  
  nomo_df <- mydf %>%
    dplyr::mutate(nomoscore = predict(a,mydf))
  
  times <- c(365,1095,1825)
  
  for (i in 1:length(times)) {
    f <- psm(Surv(OStime,OS) ~ nomoscore, data =  nomo_df, x=T, y=T, dist='lognormal') 
    ## construction of calibration curve
    cal <- calibrate(f, cmethod='KM', method="boot", u = times[i], m= 55, B = 80)
    png(filename = here::here("Plot",paste0("Nomogram_",times[i],"_year_survival.png")),res = 300,width = 6,height = 6,units = "in")
    customplot(cal,lwd=2,lty=1,
               errbar.col=ggsci::pal_nejm()(3)[1],
               xlab=paste0("Nomogram Predicted Probability of ",times[i]/365,"-Year OS"),
               ylab= paste0("Actual ",times[i]/365,"-Year OS (proportion)"),
               col=ggsci::pal_nejm()(2)[2],subtitles = FALSE)
    title(main = paste0("Calibration Curve for ",times[i]/365," Year Overall Survival"))
    dev.off()
  }
}


