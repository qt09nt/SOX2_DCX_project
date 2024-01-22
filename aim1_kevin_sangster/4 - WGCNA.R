rm(list = ls())
library(pheatmap)
library(MCMCglmm)
library(ggplot2)
library(EnhancedVolcano)
library(stringr)
library(hrbrthemes)
library(factoextra)
library(dendextend)
library(purrr)
library(gridExtra)
library(cowplot)
library(tidyverse)
library(viridis)
library(WGCNA)
library(ggsignif)
library(ggpubr)

### Run to install BioConductor packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("topGO") 
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("limma")
# BiocManager::install("WGCNA")

### installation of AnRichment package
# source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R");
# installAnRichment();

# change this to the directory you are working in
setwd("~/R projects/2020-2021/2021-06-29_AIM1_Analysis_Queenie/")

# run code containing additional functions I made
source("perseus_functions.R")
source("main_analysis_functions.R")

# module preservation and GO enrichment both take a while to run so don't need
# to run them every time, can just save results/load when necessary
RUN_MODULE_PRES = T
RUN_GO_ENRICHMENT = F

### load pre-processed data
# load FFPE data
data_jen = read.table("results/data_jen_pooled.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)
metadata_jen = read.table("results/metadata_jen_pooled.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t")

# load second FFPE dataset for further comparison
data_ugi = read.table("results/data_ugi_processed.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)
metadata_ugi = read.table("data/Ugi/metadata_brain.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t")

# load pre-processed organoid data (YFP+ samples)
data_pos = read.table("results/data_pos.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)
metadata_pos = read.table("results/metadata_pos.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t")

# specify order of groups for plotting later
metadata_pos$cell_type = factor(metadata_pos$cell_type, levels = c("SOX2", "DCX"))
metadata_jen$group = factor(metadata_jen$group, levels = c("NP", "MN"))

### Identify protein modules using WGCNA
### from: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-auto.R ###
thresh = pickSoftThreshold(t(data_jen))$powerEstimate
modules = blockwiseModules(t(data_jen), power = thresh)

# create dendrogram of protein clustering
dend = modules$dendrograms[[1]]
pdf("figures/cluster_dendrogram.pdf", width = 8, height = 4)
plotDendroAndColors(dend, modules$colors, dendroLabels = FALSE)
dev.off()

# quality control stuff
table(modules$goodGenes)
table(modules$goodSamples)

# calculate module eigengenes (summary of module expression)
eigengenes = moduleEigengenes(t(data_jen), modules$colors)
eigengene_table = eigengenes$eigengenes

## assess module preservation going from FFPE -> organoids
if(RUN_MODULE_PRES){
  
  # put data and module names in correct format
  multi = multiData(reference_dataset = t(as.matrix(data_jen)),
                    test_dataset = t(as.matrix(data_pos)))
  color = list(reference_dataset = modules$colors)
  
  # calculate module preservation
  mp = modulePreservation(multi, color, verbose = 3,
                          nPermutations = 100)
  
  # get summary preservation z-score 
  stats_z = mp$preservation$Z$ref.reference_dataset$inColumnsAlsoPresentIn.test_dataset
  stats_z = stats_z[order(stats_z[,2], decreasing = T),c(1:2)] 
  
  # get summary preservation p-value
  stats_p = mp$preservation$log.pBonf$ref.reference_dataset$inColumnsAlsoPresentIn.test_dataset
  stats_p = stats_p[rownames(stats_z), c(1:2)]
  
  stats = cbind.data.frame(stats_z,
                           log.p.Bonfsummary.pres = stats_p$log.p.Bonfsummary.pres)
  
  write.table(stats, "results/jen_to_ours_module_pres.txt", quote = F, col.names = T, row.names = T, sep = "\t")
}else{
  stats = read.table("results/jen_to_ours_module_pres.txt", quote = "", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
}

### plot module preservation stats

# remove "gold" module -> comprised of proteins not mapped to any other module
stats = stats[rownames(stats)!="gold",]

stats$name = rownames(stats)
stats$neg.log.p.bonf = -stats$log.p.Bonfsummary.pres

ggbarplot(stats, x="name", y="Zsummary.pres", fill = "grey")+
        rotate_x_text(angle = 45) +
        geom_hline(yintercept = 2, linetype = 2) +
        xlab("Module") +
        ylab("Preservation Z-score")

ggsave("figures/module_preservation_z.pdf", width = 3.0)          

ggbarplot(stats, x="name", y="neg.log.p.bonf", fill = "grey")+
  rotate_x_text(angle = 45) +
  geom_hline(yintercept = 1.30103, linetype = 2) +
  xlab("Module") +
  ylab("Preservation P-value (-log10)")

ggsave("figures/module_preservation_p.pdf", width = 3.0)          


### get the module eigengenes for our data
overlap = intersect(rownames(data_jen), rownames(data_pos))
data_overlap = matrix(nrow = nrow(data_jen), ncol = ncol(data_pos))
data_overlap = as.data.frame(data_overlap)
rownames(data_overlap) = rownames(data_jen)
colnames(data_overlap) = colnames(data_pos)

for(gene in rownames(data_jen)){
  if(gene %in% overlap){
    data_overlap[gene,] = data_pos[gene,]
  }else{
    data_overlap[gene,] = rep(NA, ncol(data_pos))
  }
}

CO_module_eigengenes = moduleEigengenes(t(data_overlap), modules$colors)
CO_eigengene_table = CO_module_eigengenes$eigengenes


### from the preserved modules, examine which ones are significantly different in NP vs MN, SOX2 vs DCX, etc.
preservation_zscore_cutoff = 3
preserved_modules = paste("ME", rownames(stats)[stats$Zsummary.pres >= preservation_zscore_cutoff], sep = "")

# remove "gold" module -> comprised of proteins not mapped to any other module
preserved_modules = preserved_modules[preserved_modules != "MEgold"]

# loop through all preserved modules and make a bunch of plots
for(module in preserved_modules){
  print(get_boxplot(eigengene_table[,module], metadata_jen$group, module,
                    adj_mult = length(preserved_modules), x_lab = "Group"))
  ggsave(paste("figures/", module, " fetal data groups.pdf"), width = 2.5)
  
  print(get_boxplot(eigengene_table[,module], metadata_jen$Neurological.Region,
                    title = module, x_lab = "Region"))
  ggsave(paste("figures/", module, " fetal data regions.pdf"), width = 4.5)
  
  print(get_boxplot(CO_eigengene_table[,module], metadata_pos$cell_type,
                    module, adj_mult = length(preserved_modules),
                    x_lab = "Cell Type"))
  ggsave(paste("figures/", module, " CO.pdf"), width = 2.5)
  
}

### Combine preserved modules with proteins differentially expressed in DCX/SOX2
### to identify targets
DCX_SOX2_tt = read.table("results/DCX_vs_SOX2_top_table.txt", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
sig = DCX_SOX2_tt$adj.P.Val <= 0.1
table(sig)

mapped_modules = modules$colors[rownames(DCX_SOX2_tt)]
preserved_modules = rownames(stats)[stats$Zsummary.pres >= preservation_zscore_cutoff]
in_preserved_pathway = mapped_modules %in% preserved_modules
table(in_preserved_pathway)

table(sig & in_preserved_pathway)
targets = DCX_SOX2_tt[sig & in_preserved_pathway,]
targets = targets[order(targets$logFC, decreasing = T),]
targets$module = modules$colors[match(rownames(targets), names(modules$colors))]

write.table(targets, "results/targets_jen.txt", sep = "\t", col.names = T, row.names = T)


### make some plots for targets:

target_plots = function(target, show_pair_pval = F){
  print(get_boxplot(data_pos[target,], metadata_pos$cell_type, title = target,
                    show_pair_pval = show_pair_pval))
  ggsave(paste("figures/", target, "_CO.pdf", sep = ""), width = 2.5)
  
  if(target %in% rownames(data_jen)){
    print(get_boxplot(data_jen[target,], metadata_jen$group, title = target,
                      show_pair_pval = show_pair_pval))
    ggsave(paste("figures/", target, "_jen_group.pdf", sep = ""), width = 2.5)
    print(get_boxplot(data_jen[target,], metadata_jen$Neurological.Region,
                      title = target, show_pair_pval = show_pair_pval))
    ggsave(paste("figures/", target, "_jen_region.pdf", sep = ""), width = 2.5)
  }
  if(target %in% rownames(data_ugi)){
    g = get_boxplot(as.numeric(data_ugi[target,]), metadata_ugi$Sample.Type,
                    title = target, show_pair_pval = show_pair_pval)
    print(g)
    ggsave(paste("figures/", target, "_UGI.pdf", sep = ""))
    
  }
}

DCX_SOX2_sig = rownames(DCX_SOX2_tt)[DCX_SOX2_tt$adj.P.Val <= 0.1]
DCX_SOX2_sig_strict = rownames(DCX_SOX2_tt)[DCX_SOX2_tt$adj.P.Val <= 0.05]

target_plots("MARCKS")
modules$colors["MARCKS"]
"MARCKS" %in% DCX_SOX2_sig_strict
# not sig in MN vs NP !!!

target_plots("RUVBL2")
modules$colors["RUVBL2"]
"RUVBL2" %in% DCX_SOX2_sig
# sig in MN vs NP

target_plots("DCLK1")
modules$colors["DCLK1"]
"DCLK1" %in% DCX_SOX2_sig_strict
# sig in MN vs NP


## Perform enrichment analysis to determine what these modules are
if(RUN_GO_ENRICHMENT){
  library(anRichment)
  entrez = convert2entrez(organism = "human", symbol = rownames(data_jen))
  GO_collection = buildGOcollection(organism = "human")
  module_color = modules$colors
  
  GO_enrichment = enrichmentAnalysis(
    classLabels = module_color, identifiers = entrez,
    refCollection = GO_collection,
    useBackground = "given",
    threshold = 1e-4,
    thresholdType = "Bonferroni",
    getOverlapEntrez = TRUE,
    getOverlapSymbols = TRUE,
    ignoreLabels = "grey");
  
  enrichment_table = GO_enrichment$enrichmentTable
  #enrichment_table = enrichment_table[enrichment_table$class %in% preserved_modules,]
  write.table(enrichment_table, "results/jen_modules_GO_enrichment.txt", sep = "\t", row.names = T, col.names = T)
  
  # need to do this because anRichment causes a conflict with WGCNA
  detach("package:anRichment", unload=TRUE)
}else{
  enrichment_table = read.table("results/jen_modules_GO_enrichment.txt", sep = "\t",
                                header = T)
}

### plot results of enrichment analysis
enrichment_table$neg.log.p = -log10(enrichment_table$pValue)

# select only the 3 preserved modules
enrichment_table = enrichment_table[enrichment_table$class %in% c("turquoise", "blue", "red"),]

enrichment_table$dataSetName = factor(enrichment_table$dataSetName,
                                      levels = enrichment_table$dataSetName)

# specify the order of the modules
enrichment_table$class = factor(enrichment_table$class,
                                levels = c("blue", "red", "turquoise"))

g = ggplot(enrichment_table, aes(x=class, y=dataSetName,
                                 size = enrichmentRatio, color = neg.log.p)) +
  #this ensures that GGplot makes a scatter plot
  geom_point(alpha=0.7) +
  #this changes the sizes of the points and creates a name for the size legend
  scale_size(range = c(2, 5), name="enrichment ratio") +
  #this assigns a particular colour pallette to be used
  scale_color_viridis(option = "C",discrete=FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Module") +
  ylab("GO Term")

g$labels$colour = "-log(p)"
g

ggsave("figures/module_GO_enrichment.pdf", width = 7, height = 6)



