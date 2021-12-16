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


dir.results <- ("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/")

#calculate differential expression between Immature and Phenol1:
label.g1 <-"Immature"
label.g2 <-"Pheno1"

immature_pheno1<- cbind(immature, pheno1)
pheno1_pheno2<-cbind(pheno1, pheno2)
immature_pheno2 <- cbind(immature, pheno2)


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
volcano_ggplot(immature_pheno2_top_table, 0.05, "Immature versus Pheno2")
volcano_ggplot(pheno1_pheno2_top_table, 0.05, "Pheno1 versus Pheno2")


########change p-value cut-off to 0.10 and keep all proteins (even insignificant ones) to get 
#full list of proteins for volcano plots
#kept zeros from original intensities matrix
#log2 normalization
#ANOVA
#fold change
#p.val adjusted with Bonferroni & also with FDR

dir.results <- ("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/ANOVA_pvalue_0.10/full list")

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

#####################################################

#Plot volcano plots

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/results/ANOVA_pvalue_0.10/full list")

immature_pheno1 <- readRDS("diff_expr_padj_all_Immature_Pheno1.rds")
immature_pheno2<-readRDS("diff_expr_padj_all_Immature_Pheno2.rds")
pheno1_pheno2<-readRDS("diff_expr_padj_all_Pheno1_Pheno2.rds")

###############plot volcano plots of significantly different proteins/genes
volcano_ggplot(immature_pheno1, 0.1, "Immature versus Pheno1")
volcano_ggplot(immature_pheno2, 0.1, "Immature versus Pheno2")
volcano_ggplot(pheno1_pheno2, 0.1, "Pheno1 versus Pheno2")


#remove vertical LOG2 FC lines - done
#remove gene labels for insignificant proteins
#check why LOG2 fold change has so many values equal to 1 ************ HIGH PRIORITY
#check the LAMB1 gene on Pheno1 vs Pheno2 - extremely low fold change values



###### remove the -Inf before doing log 2 FC calculations


#
FCGBP <- as.data.frame(immature_pheno1["FCGBP",]) 
FCGBP_t <- data.frame(t(FCGBP))   #transpose so that rows are sample names and column is protein expression value for that protein
FCGBP_t$sample <- row.names(FCGBP_t)
FCGBP_keep <- data.frame(FCGBP_t[!(FCGBP_t$FCGBP) == -Inf, ])   #filter out -Inf values

test2 <- data.frame(t(FCGBP_keep)) 
test2<-test2[1,]                   #keep just expression data 



#Re-write LOG2 FC custom function for removing -Inf values before calculating the mean/log 2 fold change for each group 
res_logfc_no_inf <- apply(immature_pheno1, 1, function(x){

  #removal of -Inf values:
  keep_values <- data.frame(x)         
  keep_values$sample <- row.names(keep_values)                 #create a column to store the sample name
  keep <- data.frame(keep_values[!(keep_values$x) == - Inf,])  #filter out -Inf values
  group <<-substr(row.names(keep), start = 1, stop = 6)
  keep2 <<-keep[,1, drop=FALSE]              #keep just the expression data & row.names(sample names)
  keep2$group <- group
  
  means <-aggregate(keep2$x, list(keep2$group), FUN=mean)   #get the mean expression, for each group
  mean.g1 <-means[1,2]
  mean.g2 <-means[2,2]
  # mean.g1=mean(x[group==unique(group)[1]])
  # mean.g2=mean(x[group==unique(group)[2]])
  
  fold.change<-(2^mean.g1-2^mean.g2) / 2^mean.g1 # expression values in log scale
  return(data.frame(fold.change=fold.change))
})

df_logfc <- data.frame(gene=rownames(immature_pheno1), bind_rows(res_logfc_no_inf))


















########################################assign a very low number to 0 values to later log transform the data:

