#Sofia organoid protein secretome data analysis based on working matrix processed using Perseus

rm(list=ls())

library(pheatmap)
library(MCMCglmm)
library(ggplot2)
library(EnhancedVolcano)
library(stringr)
library(limma)
library(factoextra)
library(dendextend)
library(purrr)
library(gridExtra)
library(cowplot)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(WGCNA)
library(RColorBrewer)
library(ggsignif)
library(colorspace)
library(ggpubr)
library(cluster)
library(factoextra)
library(DEGreport)
library(dplyr)

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr")

#load functions
source("protein_organoid_secretome_perseus_working_matrix_main_analysis_functions.R")

annotations <- read.csv("Annotations.csv")
working_data <- read.csv("Working_data.csv")

row.names(working_data) <- working_data$genes
working_data[,1]<-NULL

working_data$genes <- NULL

#get rid of the first column
annotations[,1] <- NULL

#For annotations, set sample name as row.names
row.names(annotations)<- annotations$sample 
annotations$timepoints <- as.factor(annotations$timepoints)

#match the sample names between working_data & annotation
annotations$sample <- colnames(working_data)

############# do differential expression analysis between time points - includes all cell types 

unique(annotations$timepoints)
#"Week 4"  "Week 8"  "Week 12" "Week 16" "Week 10" "Week 15" "Week 6" 

#####comparing week 4 protein expression to week 6  - (all cell types at week 4 vs all cell types at week 6)
week4_vs_week6_meta <- annotations[annotations$timepoints == "Week 4" | annotations$timepoints == "Week 6", ]

#get specific samples for week 4 & week 6 time points

week4_week6_samples <-c(week4_vs_week6_meta$sample)

#keep just the week 4 & week 6 samples from working matrix
week4_vs_week6_expr <- working_data %>% 
                            select(week4_week6_samples)
  
difference = ifelse(week4_vs_week6_meta$timepoint == "Week 4", 1, 0)
design = cbind(control = 1, difference = difference)
week_4_vs_week_6_top_table = get_top_table(week4_vs_week6_expr, design)

write.table(week_4_vs_week_6_top_table, "week_4_vs_week_6_top_table.txt", sep = "\t", quote = F, row.names = T, col.names = T)
print(volcano_plot(week_4_vs_week_6_top_table, "Organoid Protein Secretome - All Cell Types - Week 4 vs Week 6"))

########## Comparison of all cell types to all cell types (time point 1 vs timepoint 2) 
#get the differentially expressed proteins and produce a top table with logFC, pvalues, t, and adj.P.val
#then plot volcano plot

#make sure to uncomment the save top tables & print volcano plot lines within the get_timepoints function if
#running for the first time, otherwise can comment that part out if re-running
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/all cell types")

#
comparisons <- c("Week 6", "Week 8", "Week 10", "Week 12", "Week 15", "Week 16")
for (i in comparisons){
  timepoint1="Week 4"
  timepoint2=i
  tp1_vs_tp_2_top_table = get_timepoints(annotations, working_data, timepoint1= "Week 4", timepoint2 = i)
  
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("all_cell_types", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}

comparisons <- c("Week 8", "Week 10", "Week 12", "Week 15", "Week 16")
for (i in comparisons){
  get_timepoints(annotations, working_data, timepoint1 = "Week 6", timepoint2 = i)
  timepoint1="Week 6"
  timepoint2=i
  tp1_vs_tp_2_top_table = get_timepoints(annotations, working_data, timepoint1= "Week 6", timepoint2 = i)
  
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("all_cell_types", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}

comparisons <- c("Week 10", "Week 12", "Week 15", "Week 16")
for (i in comparisons){
  timepoint1="Week 8"
  timepoint2=i
  get_timepoints(annotations, working_data, timepoint1 = "Week 8", timepoint2 = i)
  
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("all_cell_types", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}

comparisons <- c("Week 12", "Week 15", "Week 16")
for (i in comparisons){
  timepoint1="Week 10"
  timepoint2=i
  get_timepoints(annotations, working_data, timepoint1 = "Week 10", timepoint2 = i)
  
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("all_cell_types", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}


####NOTE: Comparison of Week 16 to everything else: 
#Note Week 16 will show up on the right hand side of the volcano plot because it is timepoint 1
#similarly for the p ranking tables, the positive values will be for Week 16, and negative ranked values 
#will be for the earlier time points
comparisons <- c("Week 4", "Week 6", "Week 8", "Week 10", "Week 12", "Week 15")
for (i in comparisons){
  timepoint1="Week 16"
  timepoint2=i
  tp1_vs_tp_2_top_table = get_timepoints(annotations, working_data, timepoint1 = "Week 16", timepoint2 = i)
  
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("all_cell_types", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}


############################

#try with just one cell line:

#keep just SOX2
SOX2_meta <- annotations[annotations$cell_type == "SOX2",]

#get expression data for just SOX2 samples:
SOX2_samples <- c(SOX2_meta$sample)
SOX2_expr <- working_data[,SOX2_samples]

#plot volcano plots, get top table and also get p rank table for GSEA
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/SOX2_only")

comparisons <- c("Week 6","Week 8", "Week 12", "Week 16")
for (i in comparisons){
  timepoint1 = "Week 4"
  timepoint2 = i
  tp1_vs_tp_2_top_table = get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                       timepoint1 = timepoint1,  timepoint2 = i, cell_type = "SOX2")
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("SOX2", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  
}


comparisons <- c("Week 8", "Week 12", "Week 16")
for (i in comparisons){
  timepoint1 = "Week 6"
  timepoint2 = i
  tp1_vs_tp_2_top_table = get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 6",  timepoint2 = i, cell_type = "SOX2")
  
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("SOX2", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}

comparisons <- c("Week 12", "Week 16")
for (i in comparisons){
  timepoint1 = "Week 8"
  timepoint2 = i
  tp1_vs_tp_2_top_table = get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 8",  timepoint2 = i, cell_type = "SOX2")
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("SOX2", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}


############################# get expression data for just DCX cell types

DCX_meta <- annotations[annotations$cell_type == "DCX", ]

#get expression data for just DCX samples:
DCX_samples <- c(DCX_meta$sample)
DCX_expr <- working_data[,DCX_samples]

#see what unique time points are available:
unique(DCX_meta$timepoints)
# Week 4  Week 6  Week 8  Week 10 Week 12 Week 15

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/DCX_only")

#comparison of DCX week 4 to other individual time points
comparisons <- c("Week 6", "Week 8", "Week 10", "Week 12", "Week 15")
for (i in comparisons){
  timepoint1= "Week 4"
  timepoint2 = i
  get_timepoints_individual_cell_types(annotations = DCX_meta, working_data = DCX_expr,
                                     timepoint1= "Week 4", timepoint2=i, cell_type="DCX")
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("DCX", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}

#comparison of DCX week 6 to other individual DCX  time points 
comparisons <- c("Week 8", "Week 10", "Week 12", "Week 15")
for (i in comparisons){
  timepoint1= "Week 6"
  timepoint2 = i
  get_timepoints_individual_cell_types(annotations = DCX_meta, working_data = DCX_expr,
                                     timepoint1= "Week 6", timepoint2=i, cell_type="DCX")
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("DCX", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}

#comparison of DCX week 8 to other individual DCX  time points 
comparisons <- c("Week 10", "Week 12", "Week 15")
comparisons
for (i in comparisons){
  timepoint1 = "Week 8"
  timepoint2 = i
  get_timepoints_individual_cell_types(annotations = DCX_meta, working_data = DCX_expr,
                                       timepoint1= "Week 8", timepoint2=i, cell_type="DCX")
  
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("DCX", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}

comparisons <- c("Week 12", "Week 15")
for (i in comparisons){
  timepoint1= "Week 10"
  timepoint2 = i
  tp1_vs_tp_2_top_table = get_timepoints_individual_cell_types(annotations = DCX_meta, working_data = DCX_expr,
                                       timepoint1= "Week 10", timepoint2=i, cell_type="DCX")
  
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("DCX", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}

get_timepoints_individual_cell_types(annotations = DCX_meta, working_data = DCX_expr,
                                     timepoint1= "Week 12", timepoint2="Week 15", cell_type="DCX")


######################## comparison of timepoints for SYN cell type only:

#keep just SYN cells from annotation:
SYN_meta <- annotations[annotations$cell_type == "SYN",]

#SYN expression data
SYN_samples <- c(SYN_meta$sample)    #get sample names for SYN cell types
SYN_expr <- working_data[,SYN_samples]

#comparison of timepoints for SYN cell type:
unique(SYN_meta$timepoints)
#Week 8  Week 10 Week 12 Week 15

#comparison of week 8 to all other time points:
comparisons <- c("Week 10", "Week 12", "Week 15")

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/SYN_only")
for (i in comparisons){
  timepoint1 = "Week 8"
  timepoint2 = i
  tp1_vs_tp_2_top_table = get_timepoints_individual_cell_types(annotations = SYN_meta, working_data = SYN_expr, 
                                       timepoint1 = "Week 8", timepoint2=i, cell_type="SYN")
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("SYN", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)  

  
  }

#comparison of week 10 SYN to all other time points:
comparisons <- c("Week 12", "Week 15")
for (i in comparisons){
  timepoint1 = "Week 10"
  timepoint2 = i
  get_timepoints_individual_cell_types(annotations = SYN_meta, working_data = SYN_expr, 
                                       timepoint1 = "Week 10", timepoint2=i, cell_type="SYN")
  ####### make pscore and FC tables for GSEA
  pscore_table = get_pscore_table(tp1_vs_tp_2_top_table)
  write.table(pscore_table, paste("SYN", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
              sep = "\t", col.names = F, row.names = F, quote = F)  
  
}

#comparison of SYN week 12 to SYN week 15
timepoint1 = "Week 12"
timepoint2 = "Week 15"
SYN_week12_v_SYN_week15_top_table = get_timepoints_individual_cell_types(annotations = SYN_meta, working_data = SYN_expr, 
                                     timepoint1 = "Week 12", timepoint2="Week 15", cell_type="SYN")

pscore_table = get_pscore_table(SYN_week12_v_SYN_week15_top_table)
write.table(pscore_table, paste("SYN", timepoint1, timepoint2, "pscore.rnk", sep = "_"),
            sep = "\t", col.names = F, row.names = F, quote = F)  


