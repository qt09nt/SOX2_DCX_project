rm(list = ls())
library(VennDiagram)
library(pheatmap)
library(limma)
library(RColorBrewer)
library(ggplot2)
library(hrbrthemes)
library(factoextra)
library(Rtsne)
library(RColorBrewer)
library(colorspace)
library(ggpubr)

### comparison with UGI and Jen's data

# change this to the directory you are working in
setwd("~/R projects/2020-2021/2021-06-29_AIM1_Analysis_Queenie/")

# run code containing additional functions I made
source("perseus_functions.R")
source("main_analysis_functions.R")

### load data
# load pre-processed organoid data (YFP+ samples)
data_pos = read.table("results/data_pos.txt", stringsAsFactors = FALSE,
                      quote = "", header = TRUE, sep = "\t", check.names = FALSE)
metadata_pos = read.table("results/metadata_pos.txt", stringsAsFactors = FALSE,
                          quote = "", header = TRUE, sep = "\t")

# load pre-processed organoid data (all samples)
CO_data = read.table("results/data_processed.txt", sep = "\t", quote = "")
CO_metadata = read.table("results/metadata.txt", sep = "\t", quote = "",
                         header = T)

# load FFPE data for comparison
data_jen = read.table("results/data_jen_pooled.txt", stringsAsFactors = FALSE,
                      quote = "", header = TRUE, sep = "\t", check.names = FALSE)
metadata_jen = read.table("results/metadata_jen_pooled.txt",
                          stringsAsFactors = FALSE, quote = "", header = TRUE,
                          sep = "\t", check.names = FALSE)

### overlap of identified proteins
venn_diagram(list(rownames(data_jen), rownames(data_pos)),
             c("Fetal Brain" , "YFP+"), "JEN_YFP_overlap")

### spearman rank correlation of cell types vs regions ***(TODO: refactor)
# get the common proteins
intersection_jen = intersect(rownames(data_jen), rownames(data_pos))

# visualize the correlations between all samples
sample_cor = cor(data_pos[intersection_jen,],
                 data_jen[intersection_jen,],
                 method="spearman")
pheatmap(sample_cor)

# average together the samples from each region in the FFPE data
averages = cbind.data.frame(
  DSVZ = apply(sample_cor[,metadata_jen$`Neurological Region` == "DSVZ"], 1, mean),
  LGE = apply(sample_cor[,metadata_jen$`Neurological Region` == "LGE"], 1, mean),
  MGE = apply(sample_cor[,metadata_jen$`Neurological Region` == "MGE"], 1, mean),
  CX = apply(sample_cor[,metadata_jen$`Neurological Region` == "CX"], 1, mean),
  STRIA = apply(sample_cor[,metadata_jen$`Neurological Region` == "STRIA"], 1, mean)
)

# average together the DCX+ and SOX2+ samples
averages = cbind.data.frame(
  SOX2 = apply(averages[metadata_pos$cell_type == "SOX2",], 2, mean),
  DCX = apply(averages[metadata_pos$cell_type == "DCX",], 2, mean)
)

# visualize correlations (averaged)
pheatmap(averages, 
         color = sequential_hcl(10, h1=238, c1=10, l1=90,
                                h2=238, c2=43, l2=52,
                                rev = F, alpha = 1.0, power = 1),
         #filename = "figures/correlation_regions.pdf",
         width = 3,
         height = 5,
         cluster_cols = F,
         cluster_rows = T
         )

### Combined PCA
merged = cbind.data.frame(data_pos[intersection_jen,],
                          data_jen[intersection_jen,])

data_prcomp = prcomp(t(merged), scale = TRUE)$x
type = c(metadata_pos$cell_type, metadata_jen$`Neurological Region`)
source = c(rep("CO", ncol(data_pos)), rep("fetal", ncol(data_jen)))

pc_data = data.frame(x=data_prcomp[,2], y=data_prcomp[,3],  # PC1 separates based on dataset
                     type=type,
                     source=source)

### k-means clustering
kmeans_data = pc_data[,c("x", "y")]
clusters = kmeans(kmeans_data, 2, nstart = 25)
clusters$cluster[order(clusters$cluster)]

pdf("figures/k_means.pdf", width = 7, height = 7)
fviz_cluster(clusters, data = kmeans_data,
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_ipsum(),
             show.clust.cent = F)
dev.off()

### plot the PCA
pc_data$type = factor(pc_data$type, levels = c("DSVZ", "LGE", "MGE", "SOX2",
                                                   "DCX", "CX", "STRIA"))
labels = c(DSVZ="DSVZ",
           LGE="LGE",
           MGE="MGE",
           SOX2="SOX2+ (Organoid)",
           DCX="DCX+ (Organoid)",
           CX="CX",
           STRIA="STRIA")

# plot with custom colours and lebels
ggplot(pc_data, aes(x=x, y=y, color=type, shape=type)) + 
  geom_point(size=6) +
  scale_color_manual("Sample", 
                     values = c(DSVZ="#4E80A5",
                                LGE="#BD7691",
                                MGE="#AB92C5",
                                SOX2="#9DB4DD", 
                                DCX="#A4D49C",
                                CX="#FDC98F",
                                STRIA="#F69679"),
                     labels = labels) +                          
  scale_shape_manual("Sample", 
                     values = c(DSVZ=17,
                                LGE=17,
                                MGE=17,
                                SOX2=16, 
                                DCX=16,
                                CX=17,
                                STRIA=17),
                     labels = labels)+
  theme_bw() +
  xlab("PC2") +
  ylab("PC3") +
  theme(text = element_text(size = 20))


ggsave("figures/combined_PCA_pc_2_3.pdf", width = 7.5, height = 5)

