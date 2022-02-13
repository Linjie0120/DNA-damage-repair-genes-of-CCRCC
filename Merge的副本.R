library(tidyverse)
library(cowplot)


fig_LASSO <- ggdraw(xlim = c(0,15),ylim = c(0,5)) +
  draw_image(image = here::here("ccRCC_Project","Plot","LASSO_FIT.png"),
             x = 0,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","cvLASSO.png"),
             x = 5,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","coef_plot.tiff"),
             x = 10,y =0,width = 5,height = 5) + 
  draw_plot_label(label = LETTERS[1:3],x = c(0,5,10),y = c(5,5,5))

ggsave(plot = fig_LASSO,
       filename = here::here("ccRCC_Project","rs","LASSO.tiff"),
       width = 15,height = 5,
       device = "tiff",dpi = 600,compression = "lzw")

#Figure 2
fig2 <- ggdraw(xlim = c(0,15),ylim = c(0,11)) +
  draw_image(image = here::here("ccRCC_Project","Plot","Training_Risk.tiff"),
             x = 0,y =7,width = 5,height = 4) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Training_Survival.tiff"),
             x = 0,y =3,width = 5,height = 4) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Training_heatmap.tiff"),
             x = 0,y =0,width = 5,height = 3) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Testing_Risk.tiff"),
             x = 5,y =7,width = 5,height = 4) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Testing_Survival.tiff"),
             x = 5,y =3,width = 5,height = 4) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Testing_heatmap.tiff"),
             x = 5,y =0,width = 5,height = 3) +
  draw_image(image = here::here("ccRCC_Project","Plot","Entire_Risk.tiff"),
             x = 10,y =7,width = 5,height = 4) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Entire_Survival.tiff"),
             x = 10,y =3,width = 5,height = 4) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Entire_heatmap.tiff"),
             x = 10,y =0,width = 5,height = 3) +
  draw_plot_label(label = LETTERS[1:9],x = c(0,0,0,5,5,5,10,10,10),y = c(11,7,3,11,7,3,11,7,3))

ggsave(plot = fig2,
       filename = here::here("ccRCC_Project","rs","Risk_survival.tiff"),
       width = 15,height = 11,
       device = "tiff",dpi = 600,compression = "lzw")


fig_GSVA <- ggdraw(xlim = c(0,25),ylim = c(0,5)) +
  draw_image(image = here::here("ccRCC_Project","Plot","GSVA_Score_BASE_EXCISION_REPAIR.tiff"),
             x = 0,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","GSVA_Score_HOMOLOGOUS_RECOMBINATION.tiff"),
             x = 5,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","GSVA_Score_MISMATCH_REPAIR.tiff"),
             x = 10,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","GSVA_Score_NON_HOMOLOGOUS_END_JOINING.tiff"),
             x = 15,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","GSVA_Score_NUCLEOTIDE_EXCISION_REPAIR.tiff"),
             x = 20,y =0,width = 5,height = 5) + 
  draw_plot_label(label = LETTERS[1:5],x = c(0,5,10,15,20),y = c(5,5,5,5,5))

ggsave(plot = fig_GSVA,
       filename = here::here("ccRCC_Project","rs","GSVA.tiff"),
       width = 25,height = 5,
       device = "tiff",dpi = 600,compression = "lzw")


fig_survival <- ggdraw(xlim = c(0,10),ylim = c(0,15)) +
  draw_image(image = here::here("ccRCC_Project","Plot","Training_OS.tiff"),
             x = 0,y =10,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Training_PFS.tiff"),
             x = 5,y =10,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Testing_OS.tiff"),
             x = 0,y =5,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Testing_PFS.tiff"),
             x = 5,y =5,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Entire_OS.tiff"),
             x = 0,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Entire_PFS.tiff"),
             x = 5,y =0,width = 5,height = 5) + 
  draw_plot_label(label = LETTERS[1:6],x = c(0,5,0,5,0,5),y = c(15,15,10,10,5,5))
ggsave(plot = fig_survival,
       filename = here::here("ccRCC_Project","rs","risk_prognosis.tiff"),
       width = 10,height = 15,
       device = "tiff",dpi = 600,compression = "lzw")


fig5 <- ggdraw(xlim = c(0,20),ylim = c(0,20)) +
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_agehigh.tiff"),
             x = 0,y =15,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_agelow.tiff"),
             x = 5,y =15,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_female.tiff"),
             x = 10,y =15,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_male.tiff"),
             x = 15,y =15,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_G1-G2.tiff"),
             x = 0,y =10,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_G3-G4.tiff"),
             x = 5,y =10,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_I-II.tiff"),
             x = 10,y =10,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_III-IV.tiff"),
             x = 15,y =10,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_M0.tiff"),
             x = 0,y =5,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_M1.tiff"),
             x = 5,y =5,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_N0.tiff"),
             x = 10,y =5,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_N1.tiff"),
             x = 15,y =5,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_T1-T2.tiff"),
             x = 5,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","subgroup_survival_analysis_T3-T4.tiff"),
             x = 10,y =0,width = 5,height = 5) + 
  draw_plot_label(label = LETTERS[1:14],x = c(0,5,10,15,0,5,10,15,0,5,10,15,5,10),y = c(20,20,20,20,15,15,15,15,10,10,10,10,5,5))

ggsave(plot = fig5,
       filename = here::here("ccRCC_Project","rs","Subgroup_Survival.tiff"),
       width = 20,height = 20,
       device = "tiff",dpi = 600,compression = "lzw")



fig_mutation <- ggdraw(xlim = c(0,20),ylim = c(0,18)) +
  draw_image(image = here::here("ccRCC_Project","Plot","Low_Risk_onco_EMT.png"),
             x = 0,y =12,width = 10,height = 6) + 
  draw_image(image = here::here("ccRCC_Project","Plot","High_Risk_onco_EMT.png"),
             x = 10,y =12,width = 10,height = 6) + 
  draw_image(image = here::here("ccRCC_Project","Plot","mafcompare_onco.png"),
             x = 0,y =6,width = 8,height = 6) + 
  draw_image(image = here::here("ccRCC_Project","Plot","SETD2_Mut.png"),
             x = 8,y =6,width = 6,height = 6) + 
  draw_image(image = here::here("ccRCC_Project","Plot","PRKDC_Mut.png"),
             x = 14,y =6,width = 6,height = 6) + 
  draw_image(image = here::here("ccRCC_Project","Plot","tmb_twogroup.tiff"),
             x = 1,y =0,width = 6,height = 6) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Immune_Threpay_Percent.tiff"),
             x = 7,y =0,width = 6,height = 6) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Immune_Threpay_response.tiff"),
             x = 13,y =0,width = 6,height = 6) +  
  draw_plot_label(label = LETTERS[1:8],x = c(0,10,0,8,14,1,7,13),y = c(18,18,12,12,12,6,6,6))
ggsave(plot = fig_mutation,
       filename = here::here("ccRCC_Project","rs","Immune_response.tiff"),
       width = 20,height = 18,
       device = "tiff",dpi = 600,compression = "lzw")

figGO <- ggdraw(xlim = c(0,16),ylim = c(0,24)) +
  draw_image(image = here::here("ccRCC_Project","Plot","GO_plot.tiff"),
             x = 0,y =16,width = 8,height = 8) + 
  draw_image(image = here::here("ccRCC_Project","Plot","KEGG_plot.tiff"),
             x = 8,y =16,width = 8,height = 8) + 
  draw_image(image = here::here("ccRCC_Project","Plot","GO_Chord.png"),
             x = 0,y =8,width = 8,height = 8) + 
  draw_image(image = here::here("ccRCC_Project","Plot","KEGG_Chord.png"),
             x = 8,y =8,width = 8,height = 8) + 
  draw_image(image = here::here("ccRCC_Project","Plot","KEGG_GSVA_plot.tiff"),
             x = 0,y =0,width = 16,height = 8) + 
  draw_plot_label(label = LETTERS[1:5],x = c(0,8,0,8,0),y = c(24,24,16,16,8))
ggsave(plot = figGO,
       filename = here::here("ccRCC_Project","rs","GO.tiff"),
       width = 16,height = 24,
       device = "tiff",dpi = 600,compression = "lzw")


figNomogram <- ggdraw(xlim = c(0,15),ylim = c(0,22.5)) +
  draw_image(image = here::here("ccRCC_Project","Plot","univariate_rs.png"),
             x = 0,y =15,width = 7.5,height = 7.5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","multivariate_rs.png"),
             x = 7.5,y =15,width = 7.5,height = 7.5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Nomogram.png"),
             x = 0,y =5,width = 15,height = 10) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Nomogram_365_year_survival.png"),
             x = 0,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Nomogram_1095_year_survival.png"),
             x = 5,y =0,width = 5,height = 5) +  
  draw_image(image = here::here("ccRCC_Project","Plot","Nomogram_1825_year_survival.png"),
             x = 10,y =0,width = 5,height = 5) + 
  draw_plot_label(label = LETTERS[1:6],x = c(0,7.5,0,0,5,10),y = c(22.5,22.5,15,5,5,5))
ggsave(plot = figNomogram,
       filename = here::here("ccRCC_Project","rs","Nomogram.tiff"),
       width = 15,height = 22.5,
       device = "tiff",dpi = 600,compression = "lzw")

#Figure immune
fig2 <- ggdraw(xlim = c(0,20),ylim = c(0,10)) +
  draw_image(image = here::here("ccRCC_Project","Plot","Immune_cell_infiltration_RiskScore.tiff"),
             x = 0,y =5,width = 15,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Correlation_ImmuneMarker_Genes.tiff"),
             x = 15,y = 5,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Correlation_Score_CTLA4.tiff"),
             x = 0,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Correlation_Score_PD1.tiff"),
             x = 5,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Correlation_Score_PDL1.tiff"),
             x = 10,y =0,width = 5,height = 5) + 
  draw_image(image = here::here("ccRCC_Project","Plot","Correlation_Score_PDL2.tiff"),
             x = 15,y =0,width = 5,height = 5) +
  draw_plot_label(label = LETTERS[1:6],x = c(0,15,0,5,10,15),y = c(10,10,5,5,5,5))

ggsave(plot = fig2,
       filename = here::here("ccRCC_Project","rs","Immune_remake.tiff"),
       width = 20,height = 10,
       device = "tiff",dpi = 600,compression = "lzw")
