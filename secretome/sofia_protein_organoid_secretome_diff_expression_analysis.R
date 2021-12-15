#Sofia Secretome Differential Expression Analysis

#article on Why, When and How to Adjust Your P Values?
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099145/
#doi: 10.22074/cellj.2019.5992
##Jafari M, Ansari-Pour N. Why, When and how to adjust your P values? Cell J. 2019; 20(4): 604-607. doi: 10.22074/cellj.2019.5992.

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr")

library(readxl)
library(dplyr)
library(tidyr)

#load functions for analysis:
source("protein_secretome_diff_expr_functions.R")

prot_organoid_secretome <-readRDS("prot_organoid_secretome.rds")
tsne_clusters <- read.csv("01_tsne_clusters.csv")


#get how many zeros there are for each gene across samples
rowSums(prot_organoid_secretome == 0) 
colSums(prot_organoid_secretome == 0)

#log 2 transform the data
prot_organoid_secretome_log <- log2(prot_organoid_secretome)
#log transformation of this data leaves -Inf values where there were originally 0s

saveRDS(prot_organoid_secretome_log, "prot_organoid_secretome_log.rds")

#Run this on subsequent runs of this code:
prot_organoid_secretome_log<-readRDS("prot_organoid_secretome_log.rds")


########## Keeping the 0 values




#rename prot_organoid_secretome columns with cluster assignments:
clusters <- c(unique(tsne_clusters$cluster))


#get just the columns which belong to specific clusters and put into separate dataframes
immature <- subset_cluster(tsne_clusters, "Immature")
pheno1 <- subset_cluster(tsne_clusters, "Pheno1")
pheno2 <- subset_cluster(tsne_clusters,"Pheno2")

#recombine all clusters together
two_clusters <- cbind(immature, pheno1)
prot <- cbind(two_clusters, pheno2)   #protein expression matrix with cluster assignment in colnames

#remove two_clusters df since redundant 
rm(two_clusters)


dir.results <- ("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results")


#calculate differential expression between Immature and Phenol1:
label.g1 <-"Immature"
label.g2 <-"Pheno1"

calcDiffExpr(label.g1, label.g2)


#comparison of Phenol1 to Phenol2:
label.g1 = "Pheno1"
label.g2 = "Pheno2"

calcDiffExpr(label.g1, label.g2)

#comparison of Immature to Phenol2:
label.g1 = "Immature"
label.g2 = "Pheno2"

calcDiffExpr(label.g1, label.g2)



#Correct for multiple testing - calculate p adjusted value & 
#plot volcano plots for this before doing Enrichment analysis with Enrichr
#https://rcompanion.org/rcompanion/f_01.html


#############################################################

#read in differential expression top tables:

immature_pheno1_top_table <-readRDS("results/diff_expr_padj_Immature_Pheno1.rds")

immature_pheno2_top_table <- readRDS("results/diff_expr_padj_Immature_Pheno2.rds")

pheno1_pheno2_top_table <- readRDS("results/diff_expr_padj_Pheno1_Pheno2.rds")

library(ggplot2)
library(ggrepel)

###############plot volcano plots of significantly different proteins/genes
volcano_ggplot(immature_pheno1_top_table, 0.05, "Immature versus Pheno1")
volcano_ggplot(immature_pheno2_top_table, 0.05, "Immature versus Pheno1")
volcano_ggplot(pheno1_pheno2_top_table, 0.05, "Immature versus Pheno1")


###change p-value cut-off to 0.10 and keep all proteins (even insignificant ones)
dir.results <- ("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/ANOVA_pvalue_0.10")

#calculate differential expression between Immature and Phenol1:
label.g1 <-"Immature"
label.g2 <-"Pheno1"

calcDiffExpr2(label.g1, label.g2)


#comparison of Phenol1 to Phenol2:
label.g1 = "Pheno1"
label.g2 = "Pheno2"

calcDiffExpr2(label.g1, label.g2)

#comparison of Immature to Phenol2:
label.g1 = "Immature"
label.g2 = "Pheno2"

calcDiffExpr2(label.g1, label.g2)










########################################assign a very low number to 0 values to later log transform the data:

