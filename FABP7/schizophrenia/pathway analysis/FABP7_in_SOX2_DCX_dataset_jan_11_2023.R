#analysis of FABP7 expression in the SOX2 and DCX in Sofia's human cerebral organoids dataset

#January 11, 2023

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/")

source("R/perseus_functions.R")
source("R/main_analysis_functions.R")
source("R/protein_organoid_working_matrix_main_analysis_functions.R")

library(tidyverse)
library(viridis)
library(forcats)    #used for handling factors for gsea pathways dataframe & reordering by FDR.qvalue
library(dplyr)
library(stringr)  #for subsetting gsea results for separate Gene Ontology Molecular Functions (GOMF), Gene Ontology Cellular Components (GOCC)
#and GOBP (Biological Processes) plots
library(pheatmap)
library(MCMCglmm)
library(ggplot2)
library(EnhancedVolcano)
library(limma)
library(factoextra)
library(dendextend)
library(purrr)
library(gridExtra)
library(cowplot)
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

#################

#read in previously processed data (ie. where only YFP+ samples are kept)
data_pos <- read.table(file = "results/data_pos.txt", sep = "\t", header = TRUE, row.names = 1)
metadata_pos <- read.table(file = "results/metadata_pos.txt", sep = "\t", header = TRUE)

data_pos$genes <- row.names(data_pos)

#filter for just the row containing FABP7 expression values
FABP7 <- data_pos[data_pos$genes == "FABP7",]

FABP7_t<- t(FABP7)

FABP7_t <-as.data.frame(FABP7_t[1:13,])

row.names(FABP7_t)<-metadata_pos$sample

combined <-cbind(FABP7_t, metadata_pos)

colnames(combined)[1]<- "FABP7"

combined$tp <- factor(combined$tp, levels = c('1', '2', '3'))

combined$FABP7 <- as.numeric(combined$FABP7)

#boxplot of FABP7 expression in SoX2 populations and DCX populations
combined %>%
  ggplot(aes(x=tp, y=FABP7, fill=cell_type))+
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size = 0.4, alpha=0.9) +
  theme(
    plot.title = element_text(size=11)
  ) +
  ggtitle("FABP7 Expression in SOX2 and DCX organoids")

#jitter version which includes individual sample points
ggplot(combined, aes(x=tp, y=FABP7, fill=cell_type)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0.1) +
  theme(
    plot.title = element_text(size=11)
  ) + 
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines and grey background
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        
  ) +
  
  #ggtitle("RUVBL2 Expression in ganglionic eminence and cortical regions") +
  labs(title= "FABP7 Expression in SOX2 and DCX Organoids",
       y="FABP7 expression (LFQ intensity)",
       x="Time Point")
