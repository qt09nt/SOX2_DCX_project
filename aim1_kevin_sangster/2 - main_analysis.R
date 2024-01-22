rm(list = ls())
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

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("topGO") 
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("limma")
# BiocManager::install("WGCNA")

# change this to the directory you are working in
setwd("~/R projects/2020-2021/2021-06-29_AIM1_Analysis_Queenie/")

# run code containing additional functions I made
source("perseus_functions.R")
source("main_analysis_functions.R")

### load data
# load pre-processed organoid data (all samples)
CO_data = read.table("results/data_processed.txt", sep = "\t", quote = "")
CO_metadata = read.table("results/metadata.txt", sep = "\t", quote = "",
                         header = T)

# load pre-processed organoid data (YFP positive samples)
data_pos = read.table("results/data_pos.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)
metadata_pos = read.table("results/metadata_pos.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t")

### generate principle component plots (scale=TRUE isn't necessary but is advised in the function description)
data_prcomp = prcomp(t(data_pos), scale = TRUE)$x

# plot (for example PC1 by PC2)
pc_data = data.frame(x=data_prcomp[,1], y=data_prcomp[,2], tp=as.character(metadata_pos$tp),
                     cell_type=metadata_pos$cell_type)
ggplot(pc_data, aes(x=x, y=y, color=cell_type, shape=cell_type)) + 
  geom_point(size=6) +
  scale_color_manual("Cell Type",
                     values = c("#A4D49C", "#4E80A5")) +
  scale_shape("Cell Type") +
  theme_bw() +
  xlab("PC1") +
  ylab("PC2") +
  theme(text = element_text(size = 14)) 
  

ggsave("figures/PCA.pdf", width = 5, height = 4)

### do kmeans clustering
fviz_nbclust(t(data_pos), kmeans, method = "silhouette")
clusters = kmeans(t(data_pos), 2, nstart = 25)
clusters$cluster[order(clusters$cluster)]

fviz_cluster(clusters, data = t(data_pos),
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw(),
             show.clust.cent = F)

#### visualize column correlations with heatmap

# generate column correlations
data_cor = cor(data_pos, method = "pearson")

# create custom sample names for the heatmap
custom_names = paste(metadata_pos$cell_type, "+ ",
                     c("wk 4", "wk 6", "wk 8")[metadata_pos$tp], sep = "")

# annotate the samples by cell type
annotation_col = data.frame(metadata_pos$cell_type)
colnames(annotation_col) = c("Cell Type")
rownames(annotation_col) = colnames(data_cor)

# pick custom colours for the cell type annotations
annotation_colours = list()
annotation_colours[["Cell Type"]] = c(DCX="#A4D49C", SOX2="#4E80A5")


pheatmap(data_cor, cluster_rows = T, cluster_cols = T,
         labels_row = custom_names, labels_col = custom_names,
         annotation_col = annotation_col,
         #filename = "figures/sample_correlation_heatmap.pdf",
         width = 7.5, height = 6,
         annotation_colors = annotation_colours,
         fontsize = 14)

### Differential expression analysis between cell types
difference = ifelse(metadata_pos$cell_type == "DCX", 1, 0)
design = cbind(control=1, difference=difference)
DCX_SOX2_top_table = get_top_table(data_pos, design)

write.table(DCX_SOX2_top_table, "results/DCX_vs_SOX2_top_table.txt", sep = "\t", quote = F, row.names = T, col.names = T)
print(volcano_plot(DCX_SOX2_top_table, "DCX vs SOX2"))

# make pscore and FC tables for GSEA
pscore_table = get_pscore_table(DCX_SOX2_top_table)
write.table(pscore_table, "results/DCX_vs_SOX2_pscore.rnk", sep = "\t", col.names = F, row.names = F, quote = F)


### Differential expression analysis between DCX timepoints

FDR_cutoff = 0.1
sig_between_tp = c()
for(cell_type in c("DCX", "SOX2")){
  for(timepoint in c(1, 2, 3)){
    
    # since we are missing SOX2 timepoint 3, skip
    if(cell_type == "SOX2" && timepoint == 3){
      break
    }
    
    # get differentially expressed proteins
    top_table = get_tp_top_table(data_pos, metadata_pos, cell_type, timepoint)
    write.table(top_table, paste("results/", cell_type, timepoint, "top_table.txt", sep = "_"),
                sep = "\t", quote = F, row.names = T, col.names = T)
    
    # get pscore rankings for GSEA
    pscore_table = get_pscore_table(top_table)
    write.table(pscore_table, paste("results/", cell_type, timepoint, "pscore.rnk", sep = "_"),
                sep = "\t", col.names = F, row.names = F, quote = F)
    
    # visualize with volcano plot
    print(volcano_plot(top_table, paste(cell_type, "tp", timepoint)))
    
    # keep track of signficicant proteins
    sig_between_tp = c(sig_between_tp,
                       rownames(top_table)[top_table$adj.P.Val < FDR_cutoff])
  }
}

### Visualize significant proteins from differential expression analysis (DCX timepoints)

# create custom sample names for the heatmap
custom_names = paste("DCX+ ",
                     c("wk 4", "wk 6", "wk 8")[metadata_pos$tp[metadata_pos$cell_type == "DCX"]],
                     sep = "")

# annotate the samples by timepoint
data_DCX = data_pos[,metadata_pos$cell_type == "DCX"]
annotation_col = data.frame(metadata_pos$tp[metadata_pos$cell_type == "DCX"])
colnames(annotation_col) = c("Timepoint (DCX+)")
rownames(annotation_col) = colnames(data_DCX)
annotation_col[["Timepoint (DCX+)"]] = as.character(annotation_col[["Timepoint (DCX+)"]])

# pick custom colours for the timepoint annotations
annotation_colours = list()
annotation_colours[["Timepoint (DCX+)"]] = c("1"="#A4D49C",
                                             "2"="#FFF79F", 
                                             "3"="#9DB4DD")

pheatmap(data_DCX[sig_between_tp,],
         #filename = "figures/DCX_tp_heatmap.pdf",
         labels_col = custom_names,
         cluster_cols = F,
         annotation_col = annotation_col,
         annotation_colors = annotation_colours,
         width = 6,
         height = 8)

