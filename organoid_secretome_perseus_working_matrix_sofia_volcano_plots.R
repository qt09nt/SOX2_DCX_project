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

library(ggplot2)

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


##########Visualizing the GSEA results for timepoint comparisons (includes all cell types: SYN, SOX2 and DCX)

#GSEA pre-ranked was conducted on the proteome secretome timepoint comparisons p rank score files
#with databases searched in Human_AllPathways_June_01_2017_symbol.gmt, Hallmarks(ha.lla.v7.4.symbols.gmt), 
#Kegg(c2.cp.kegg.v7.4.symbols.gmt), reactome(c2.cp.reactome.v7.4.gymbols.gmt), wikipathways(c2.cp.wikipathways.v7.4.symbols.gmt),
#go bp(c5.go.bp.v7.4.symbols.gmt), go cc(c5.go.cc.v7.4.symbols.gmt), go mf(c5.go.mf.v7.4.symbols.gmt) databases

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/all cell types/GSEA/all_cell_types/all_cell_types_week_4_vs_week_6.GseaPreranked.1641334710091")

#load in the gesea_report_for_na_pos (positive enrichment for week 4 compared to week 6):
week_4_vs_week_6_na_pos_gsea <-read.table(file = 'gsea_report_for_na_pos_1641334710091.tsv', sep = '\t', header = TRUE)
colnames(week_4_vs_week_6_na_pos_gsea)[1] <-"Pathways"     #rename first column of GSEA pos results to "Pathways"

#filter the GSEA results dataframe to keep just the pathways which have a NOM p-value of less than 0.065 and a 
#FDR q values less than 0.1
filtered <- week_4_vs_week_6_na_pos_gsea[week_4_vs_week_6_na_pos_gsea$`NOM.p.val` <= 0.07 && week_4_vs_week_6_na_pos_gsea$`FDR.q.val` <= 0.09, ]
filtered$timepoint <- "week 4" 

#load in the gsea report for na_neg (negative enrichment for week 4 compared to week 6)
week_4_vs_week_6_na_neg_gsea <- read.table(file = 'gsea_report_for_na_neg_1641334710091.tsv', sep = '\t', header = TRUE)
colnames(week_4_vs_week_6_na_neg_gsea)[1] <- "Pathways"   ##rename first column of GSEA pos results to "Pathways"
filtered_timepoint_2 <- week_4_vs_week_6_na_neg_gsea[week_4_vs_week_6_na_neg_gsea$`NOM.p.val` <= 0.07 && week_4_vs_week_6_na_neg_gsea$`FDR.q.val` <= 0.09, ]
filtered_timepoint_2$timepoint <- "week 6"

#merge the GSEA positive and negative dataframes together
merged <- rbind(filtered, filtered_timepoint_2)

ggplot(merged, aes(x = timepoint, y = Pathways, size = NES, color = FDR.q.val )) +
  geom_point(alpha = 0.7) + #scatterplot
  scale_size(range = c(2, 5), name = "NES") +  #change the size of the points and create a name for the size legend
  scale_color_viridis(discrete = FALSE)

#GSEA plot of Week 4 vs Week 6 comparisons
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/all cell types/GSEA/all_cell_types_week_4_vs_week_6.GseaPreranked.1640185266242")
plot_gsea_timepoints(gsea_report_na_pos = 'gsea_report_for_na_pos_1640185266242.tsv',  
                     gsea_report_for_na_neg = 'gsea_report_for_na_neg_1640185266242.tsv',
                     timepoint1 = "week 4", timepoint2 = "week 6")

#GSEA plot of Week 4 vs Week 8 comparisons
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/all cell types/GSEA/all_cell_types/all_cell_types_week_4_vs_week_8.GseaPreranked.1641340848426")
plot_gsea_timepoints(gsea_report_na_pos = 'gsea_report_for_na_pos_1641340848426.tsv',
    gsea_report_for_na_neg = 'gsea_report_for_na_neg_1641340848426.tsv', timepoint1 = "week 4", timepoint2 = "week 8")

#week 4 vs week 10
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/all cell types/GSEA/all_cell_types/all_cell_types_week_4_vs_week_10.GseaPreranked.1641341772629")
plot_gsea_timepoints(gsea_report_na_pos = 'gsea_report_for_na_pos_1641341772629.tsv',
                     gsea_report_for_na_neg = 'gsea_report_for_na_neg_1641341772629.tsv', timepoint1 = "week 4", timepoint2 = "week 10")

#week 4 vs week 12
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/all cell types/GSEA/all_cell_types/all_cell_types_week_4_vs_week_12.GseaPreranked.1641344672042")
plot_gsea_timepoints(gsea_report_na_pos = 'gsea_report_for_na_pos_1641344672042.tsv',
                     gsea_report_for_na_neg = 'gsea_report_for_na_neg_1641344672042.tsv', timepoint1 = "week 4", timepoint2 = "week 12")

