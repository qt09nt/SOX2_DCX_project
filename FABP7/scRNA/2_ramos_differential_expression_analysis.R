#Differential expression of FABP7 higher expression clusters in Ramos et al (2022) GSM6720853
#May 2 2023


library(EnhancedVolcano)

setwd("C:/Users/Diamandis Lab II/Documents/Queenie/ramos/GSM6720853")

clusterFABP7_low_log2.markers <-readRDS("clusterFABP7_low_log2.markers.RDS")

##############Alternate Ways to conduct Differential Expression Analysis
##### using DESeq2
# install.packages(airway)
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("airway")
# 
# library(airway)
# library(magrittr)
# 
# data('airway')
# airway$dex %<>% relevel('untrt')
# 
# Annotate the Ensembl gene IDs to gene symbols:
#   
# ens <- rownames(airway)
# 
# library(org.Hs.eg.db)
# symbols <- mapIds(org.Hs.eg.db, keys = ens,
#                   column = c('SYMBOL'), keytype = 'ENSEMBL')
# symbols <- symbols[!is.na(symbols)]
# symbols <- symbols[match(rownames(airway), names(symbols))]
# rownames(airway) <- symbols
# keep <- !is.na(rownames(airway))
# airway <- airway[keep,]
# 
# Conduct differential expression using DESeq2 in order to create 2 sets of results:
#   
#   library('DESeq2')
# 
# dds <- DESeqDataSet(airway, design = ~ cell + dex)
# dds <- DESeq(dds, betaPrior=FALSE)
# res <- results(dds,
#                contrast = c('dex','trt','untrt'))
# res <- lfcShrink(dds,
#                  contrast = c('dex','trt','untrt'), res=res, type = 'normal')





######################## Pathway Analysis 
#https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html

# The ReactomeGSA package is a client to the web-based Reactome Analysis System.
# Essentially, it performs a gene set analysis using the latest version of the Reactome 
# pathway database as a backend.
# 
# This vignette shows how the ReactomeGSA package can be used to perform a pathway analysis
# of cell clusters in single-cell RNA-sequencing data.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require(ReactomeGSA))
  BiocManager::install("ReactomeGSA")
#> Loading required package: ReactomeGSA

# install the ReactomeGSA.data package for the example data
if (!require(ReactomeGSA.data))
  BiocManager::install("ReactomeGSA.data")

library(ReactomeGSA)
library(ReactomeGSA.data)

if (!require(ReactomeGSA.data))
  BiocManager::install("ReactomeGSA.data")
