library(tidyverse)
library(tidymodels)
library(msigdbr)
library(limma)
library(openxlsx)

rm(list = ls())
#Load ccRCC data
load(file = here::here("TCGA_GTEX_CCLEdata","KIRC","KIRC_Data.RData"))
source(here::here("Code","Function.R"))
load(here::here("TCGA_GTEX_CCLEdata","Pancancer_merge_survival.RData"))

#mismatch repair defects
#HRR，NHEJ
kegg_db <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG") %>%
  dplyr::filter(str_detect(gs_name,pattern = "REPAIR|JOINING|HOMOLOGOUS|MICROHOMOLOGY"))

write.xlsx(kegg_db,here::here("ccRCC_Project","DataFiles","DNA_Repair_Genes.xlsx"),overwrite = T)

dnaRepair <- kegg_db %>% dplyr::distinct(gene_symbol) %>% dplyr::pull(gene_symbol)
#DEG
{
  KIRC_df <- KIRC_expr %>% dplyr::select(KIRC_pheno_df$TCGA_id)
  group_list <- KIRC_pheno_df$type
  group_list <- factor(group_list,levels = c("Tumor","Normal"))
  
  #矩阵
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  head(design)
  rownames(design)=colnames(KIRC_df)
  head(design)
  identical(rownames(design),colnames(KIRC_df))
  # 比较矩阵
  contrast.matrix <- makeContrasts(Tumor-Normal,
                                   levels=design)
  fit <- lmFit(KIRC_df, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)   #该GSE数据集只有肿瘤与正常组织，可以不要contrast.matrix
  fit2 <- eBayes(fit2)
  allDiffgenes<-topTable(fit2,adjust="BH",number=Inf,coef=1) %>%
    tibble::rownames_to_column(var = "Symbol")
  
  deg_repair <- allDiffgenes %>% dplyr::filter(Symbol %in% dnaRepair) %>%
    dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 0.5)
  
  degs_list <- list(degs = allDiffgenes,deg_repair = deg_repair)
  write.xlsx(degs_list,here::here("ccRCC_Project","DataFiles","DEGs.xlsx"),overwrite = T)
}

survival_rs <- gene_list_survival_analysis(genelist = deg_repair$Symbol,type = "KIRC")
survival_rs_repair <- gene_list_survival_analysis(genelist = intersect(dnaRepair,allDiffgenes$Symbol),type = "KIRC")

survival_list <- list(degs_survival = survival_rs,repair_survival = survival_rs_repair)
write.xlsx(survival_list,here::here("ccRCC_Project","DataFiles","survival_list.xlsx"),overwrite = T)

save(survival_rs,allDiffgenes,dnaRepair,file = here::here("ccRCC_Project","RData","DNA_Repair.RData"))
