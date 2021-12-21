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


############# do differential expression analysis between time points - includes all cell types 

unique(annotations$timepoints)
#"Week 4"  "Week 8"  "Week 12" "Week 16" "Week 10" "Week 15" "Week 6" 

#match the sample names between working_data & annotation
annotations$sample <- colnames(working_data)

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

########## write this into a function:

get_timepoints(annotations, working_data, timepoint1 = "Week 4", "Week 8")
get_timepoints(annotations, working_data, timepoint1 = "Week 4", "Week 10")
get_timepoints(annotations, working_data, timepoint1 = "Week 4", "Week 12")
get_timepoints(annotations, working_data, timepoint1 = "Week 4", "Week 15")

get_timepoints(annotations, working_data, timepoint1 = "Week 6", "Week 8")
get_timepoints(annotations, working_data, timepoint1 = "Week 6", "Week 10")
get_timepoints(annotations, working_data, timepoint1 = "Week 6", "Week 12")
get_timepoints(annotations, working_data, timepoint1 = "Week 6", "Week 15")

get_timepoints(annotations, working_data, timepoint1 = "Week 8", "Week 10")
get_timepoints(annotations, working_data, timepoint1 = "Week 8", "Week 12")
get_timepoints(annotations, working_data, timepoint1 = "Week 8", "Week 15")

get_timepoints(annotations, working_data, timepoint1 = "Week 10", "Week 12")
get_timepoints(annotations, working_data, timepoint1 = "Week 10", "Week 15")

############################

#try with just one cell line:

#keep just SOX2
SOX2_meta <- annotations[annotations$cell_type == "SOX2",]

#get expression data for just SOX2 samples:
SOX2_samples <- c(SOX2_meta$sample)
SOX2_expr <- working_data[,SOX2_samples]

get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 4",  timepoint2 = "Week 6", cell_type = "SOX2")
get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 4",  timepoint2 = "Week 8", cell_type = "SOX2")
get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 4",  timepoint2 = "Week 12", cell_type = "SOX2")
get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 4",  timepoint2 = "Week 16", cell_type = "SOX2")
get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 6",  timepoint2 = "Week 8", cell_type = "SOX2")
get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 6",  timepoint2 = "Week 12", cell_type = "SOX2")
get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 6",  timepoint2 = "Week 16", cell_type = "SOX2")
get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 8",  timepoint2 = "Week 12", cell_type = "SOX2")
get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 8",  timepoint2 = "Week 16", cell_type = "SOX2")
get_timepoints_individual_cell_types(annotations = SOX2_meta, working_data = SOX2_expr, 
                                     timepoint1 = "Week 12",  timepoint2 = "Week 16", cell_type = "SOX2")


##########

