#Allen Human Brain Atlas Developmental Transcriptome RNA Seq data Analysis 

#Data downloaded from https://www.brainspan.org/static/download.html RNA-Seq Gencode v10 summarized to genes

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Allen Brain Atlas/BrainSpan/data/developmental transcriptome dataset/genes_matrix_csv_2")

library(tidyverse)
library(dplyr)
library(stringr)
library(janitor)
library(ggplot2)
library(viridis)
library(data.table)
library(ggpubr)
library(broom)
library(AICcmodavg)
library(lme4)        #mixed effects model
library(mice)
library(mitml)

#this is the gene expression data in TPM format (no need to re-run the calculations for generating the TPM values below
#if it can load the TPM data here)
TPM <-readRDS("TPM.rds")

# data <- read.table("expression_matrix.txt", header=T, row.names=1, sep ="\t")
data <- read.csv(file = 'expression_matrix.csv', header=F)  #expression values in RPKM
columns_meta <- read.csv(file = 'columns_metadata.csv', header = T)
rows_meta <- read.csv(file = 'rows_metadata.csv', header = T)

data <- data[2:525]

#RNA Seq expression data is currently in RPKM (reads per kilobase per million) which is not recommended 
#for comparison of gene expression between different samples, since the total 
#number of RPKM/FPKM normalized counts for each sample will be different. Therefore, you cannot
#compare the normalized counts for each gene equally between samples.

#source: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html


# A different RNA seq unit is TPM, which represents the number of transcripts of feature 
#found in a total number of one million transcripts

#To convert from RPKM to TPM (transcripts per kilobase million), divide each feature 
#(gene expression value in RNA seq RPKM) by the sum of the RPKM values of all features (genes) in the sample
#and multiplying by one million. 
#
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7373998/
#Zhao, S., Ye, Z., & Stanton, R. (2020). Misuse of RPKM or TPM normalization when comparing across
#samples and sequencing protocols. RNA (New York, N.Y.), 26(8), 903-909. https://doi.org/10.1261/rna.074922.120

#find the sum of each of the data columns (samples) (RPKM_sums)
RPKM_sums <- colSums(data)

#make a copy of data
TPM <- data

#function for finding sum of sample column:
RPKM_sum <- function(column){
  rpkm_total <- sum(column)
  rpkm_total <<- rpkm_total
}

#function for dividing by rpkm_total:
divide <- function (value){
  v <- value/rpkm_total 
  return(v)
}

#function for multiplying by 1 million
multiply_million <- function (x){
  million <- x*10^6
  million <<- million
}

#convert RNA Seq expression data from RPKM to TPM:
for (c in 1:ncol(TPM)){
  sample <- TPM[,c]   #extract column c from TPM dataframe
  rpkm_total <- RPKM_sum(sample)      #find total RPKM for that sample (column)
  division <- map(sample, divide)     #divide every value in the sample column by the sum of RPKM
  new <- map(division, multiply_million) #multiply every value in sample by 10^6
  new <- unlist(new)
  TPM[c] <-new   #replace converted TPM values for the sample values of that column
}

#merge gene symbols from rows_meta with TPM expression data:
TPM$gene_symbol <- rows_meta$gene_symbol

#write the RNA Seq expression data (TPM) into a csv file so that it can be read later &
#don't have to re-run RPKM to TPM conversion every time you run the code:
write.csv(TPM, 'TPM.csv')

#######################################################
#read in RNA expression file which has been converted to TPM units (this can be run in subsequent re-running of this code
#after the TPM.csv has been previously generated)
TPM <- read.csv(file = 'TPM.csv', header = T)

#filter for RUVBL2
RUVBL2 <- TPM[TPM$gene_symbol == "RUVBL2",]

#check for GFAP expression
GFAP <- TPM2[TPM2$gene_symbol == "GFAP", ]

#Oct 20, 2022 filter for FABP7, FABP3, FABP5, VRK1
FABP7 <-TPM[TPM$gene_symbol == "FABP7", ]

#keep just the columns with RNA expression data values 
RUVBL2 <- RUVBL2[,2:525]
RUVBL2 <- t(RUVBL2)

GFAP <- GFAP[,2:525]
GFAP <- t(GFAP)

FABP7 <- FABP7[,2:525]
FABP7 <- t(FABP7)

#remove the "V" from row names in data dataframe:
row.names(RUVBL2) <- row.names(RUVBL2) %>% str_replace("V","")
RUVBL2 <- as.data.frame(RUVBL2)
colnames(RUVBL2)[1] <- "RUVBL2"

#remove the "V" from row names in data dataframe:
row.names(GFAP) <- row.names(GFAP) %>% str_replace("X","")
GFAP <- as.data.frame(GFAP)
colnames(GFAP)[1] <- "GFAP"

##remove the "X" from row names in data dataframe:
row.names(FABP7) <- row.names(FABP7) %>% str_replace("X","")
FABP7 <- as.data.frame(FABP7)
colnames(FABP7)[1] <- "FABP7"

#convert the row names of RUVBL2 data frame to same row names (sample numbers) as those in column_meta data frame
#so that they would match in order to merge the 2 data frames
row.names(RUVBL2) <- row.names(columns_meta)

row.names(GFAP) <- row.names(columns_meta)

row.names(FABP7) <- row.names(columns_meta)

#combine RUVBL2 expression data with meta data
combined <- cbind(RUVBL2, columns_meta)

combined <- cbind(GFAP, columns_meta)

combined <- cbind(FABP7, columns_meta)

##############################################
#re-plot RUVBL2 for cortex and ganglionic eminences for prenatal 8pcw to 24 pcw:

#keep only cortex and ganglionic eminence samples:
RUVBL2 <- dplyr::filter(combined, grepl('ganglionic eminence|cortex', structure_name))

GFAP <- dplyr::filter(combined, grepl('ganglionic eminence|cortex', structure_name))

FABP7 <- dplyr::filter(combined, grepl('ganglionic eminence|cortex', structure_name))

#November 23 2022
#filter for FABP7 data for specific structures: occipital neocortex, orbital frontal cortex, parietal neocortex, 
#cerebellar cortex, caudal ganglionic eminences, lateral ganglionic eminences
FABP7_nov <- combined[combined$structure_acronym %in% c("Ocx", "OFC", "PCx", "CBC", "CGE", "LGE", "MGE"),]

FABP7_8_9_pcw <- FABP7_nov[FABP7_nov$age %in% c("8 pcw", "9 pcw"),]


#set age as a factor
RUVBL2$age <- as.factor(RUVBL2$age)

GFAP$age <- as.factor(GFAP$age)
FABP7$age <- as.factor(FABP7$age)
combined$age <-as.factor(combined$age)
FABP7_8_9_pcw$age <- as.factor(FABP7_8_9_pcw$age)

#filter for prenatal data 8pcw to 24 pcw:
RUVBL2_prenatal <- RUVBL2[RUVBL2$age %in% c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw"),]
GFAP_prenatal <- GFAP[GFAP$age %in% c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw"),]

GFAP_adult <- GFAP[GFAP$age %in% c("18 yrs", "19 yrs", "21 yrs", "23 yrs", "30 yrs", "36 yrs", "37 yrs", "40 yrs"),]
  


  
#specify the order of ages:
RUVBL2_prenatal$age <- ordered(RUVBL2_prenatal$age, levels = c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw"))
GFAP_prenatal$age <- ordered(GFAP_prenatal$age, levels = c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw"))

GFAP_adult$age <- ordered(GFAP_adult$age, levels = c("18 yrs", "19 yrs", "21 yrs", "23 yrs", "30 yrs", "36 yrs", "37 yrs", "40 yrs"))

FABP7_8_9_pcw$age <- ordered(FABP7_8_9_pcw$age, levels = c("8 pcw", "9 pcw"))

ggplot(data = RUVBL2_prenatal, aes(x = structure_name,  y = RUVBL2, group = structure_name, colour = structure_name)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        ) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle='RUVBL2 in ganglionic eminences and cortex at 8pcw to 24 pcw',
       x = 'structure',
       y = 'RUVBL2 expression - RNA Seq (TPM)')


ggplot(data = GFAP_prenatal, aes(x = structure_name,  y = GFAP, group = structure_name, colour = structure_name)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  ) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle='GFAP in ganglionic eminences and cortex at 8pcw to 24 pcw',
       x = 'structure',
       y = 'GFAP expression - RNA Seq (TPM)')


ggplot(data = GFAP_adult, aes(x = structure_name,  y = GFAP, group = structure_name, colour = structure_name)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  ) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle='GFAP in ganglionic eminences and cortex from 18 years to 40 years',
       x = 'structure',
       y = 'GFAP expression - RNA Seq (TPM)')

########################



#extract the RNA expression data for the other genes of interest:
genes <- TPM[TPM$gene_symbol %in% c("SOX2","GABRA1", "PAX6", "MAP2", "DCX", "SRSF9", "TCERG1", "HNRNPD", "CPSF7"),]

row.names(genes) <- genes$gene_symbol
genes <- genes[,1:524]

genes_t <- t(genes)
#remove the "V" from rownames of genes
row.names(genes_t) <- row.names(genes_t) %>% str_replace("V","")

combined <- cbind(genes_t, combined)

#filter for samples from 8 pcw and 9 pcw only:
combined_ge <- combined[combined$age %in% c("8 pcw", "9 pcw"),]

#remove the "pcw" from age values:
combined_ge$age <- gsub("pcw", "", as.character(combined_ge$age))
combined_ge$age <- as.numeric(combined_ge$age)

ge_cortical_89pcw <- combined_ge %>% 
     filter(structure_name %in% c("lateral ganglionic eminence", "medial ganglionic eminence",
                                  "caudal ganglionic eminence", "dorsolateral prefrontal cortex",
                                  "orbital frontal cortex", "parietal neocortex", "occipital neocortex"))

#for ganglionic and cortex dataframe at 8 pcw and 9 pcw, add another column that says whether it
#is ganglion or cortex
ge_cortical_89pcw$group <- ifelse(grepl("GE", ge_cortical_89pcw$structure_acronym), "ganglionic eminence", "cortex")

#facet boxplots of ganglionic eminence vs cortex at 8pcw and 9pcw
r <- ggplot(data = ge_cortical_89pcw, aes(x = group,  y = RUVBL2, group = group, colour = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle='RUVBL2 expression in ganglionic eminence and cortex at 8pcw and 9pcw',
       x = 'structure',
       y = 'RUVBL2 expression (TPM)')

#Adding p value: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
#  Add p-value
# Change method
r + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))

#ANOVA test
two.way <- aov(RUVBL2 ~ group + age, data = ge_cortical_89pcw)
summary(two.way)
# 
#             Df Sum Sq Mean Sq F value Pr(>F)
# group        1  115.1   115.1   0.406  0.537
# age          1  123.8   123.8   0.437  0.522
# Residuals   11 3116.2   283.3 

#facet boxplots of ganglionic eminence vs cortex at 8pcw and 9pcw
s <- ggplot(data = ge_cortical_89pcw, aes(x = group,  y = SOX2, group = group, colour = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle='SOX2 expression in ganglionic eminence and cortex at 8pcw and 9pcw',
       x = 'structure',
       y = 'SOX2 expression (TPM)')

#Adding p value: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
#  Add p-value
# Change method
s + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))

#GABRA1
gene <- "GABRA1"
g <- ggplot(data = ge_cortical_89pcw, aes(x = group,  y = GABRA1, group = group, colour = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #axis.text.y = element_text(face="bold", color="black")
        ) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle= paste0(toString(gene), ' expression in ganglionic eminence and cortex at 8pcw and 9pcw'),
       x = 'structure',
       y = paste0(toString(gene),' expression (TPM)'))

g + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))

gene_names <- c(colnames(ge_cortical_89pcw)[1:10])
gene <- gene_names[1]

#PAX6
p <- ggplot(data = ge_cortical_89pcw, aes(x = group,  y = PAX6, group = group, colour = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #axis.text.y = element_text(face="bold", color="black")
  ) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle= paste0(toString(gene), ' expression in ganglionic eminence and cortex at 8pcw and 9pcw'),
       x = 'structure',
       y = paste0(toString(gene),' expression (TPM)'))

p + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))

#DCX
gene <- gene_names[3]
d <- ggplot(data = ge_cortical_89pcw, aes(x = group,  y = DCX, group = group, colour = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #axis.text.y = element_text(face="bold", color="black")
  ) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle= paste0(toString(gene), ' expression in ganglionic eminence and cortex at 8pcw and 9pcw'),
       x = 'structure',
       y = paste0(toString(gene),' expression (TPM)'))

d + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))

#MAP2
gene <- gene_names[4]

m <- ggplot(data = ge_cortical_89pcw, aes(x = group,  y = MAP2, group = group, colour = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #axis.text.y = element_text(face="bold", color="black")
  ) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle= paste0(toString(gene), ' expression in ganglionic eminence and cortex at 8pcw and 9pcw'),
       x = 'structure',
       y = paste0(toString(gene),' expression (TPM)'))

m + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))

#SRSF9
gene <- gene_names[5]

sr <- ggplot(data = ge_cortical_89pcw, aes(x = group,  y = MAP2, group = group, colour = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #axis.text.y = element_text(face="bold", color="black")
  ) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle= paste0(toString(gene), ' expression in ganglionic eminence and cortex at 8pcw and 9pcw'),
       x = 'structure',
       y = paste0(toString(gene),' expression (TPM)'))

sr + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))

#TCERG1
gene <- gene_names[6]


  t <- ggplot(data = ge_cortical_89pcw, aes(x = group,  y = TCERG1, group = group, colour = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
    facet_wrap(~ age) +  #separate the boxplots by age
    theme_bw() + #set the background to white
    theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          #axis.text.y = element_text(face="bold", color="black")
    ) +
    labs(title = 'BrainSpan Transcriptome Dataset',
         subtitle= paste0(toString(gene), ' expression in ganglionic eminence and cortex at 8pcw and 9pcw'),
         x = 'structure',
         y = paste0(toString(gene),' expression (TPM)'))
  
  t + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))

#HNRNPD
gene <- gene_names[7]  

h <- ggplot(data = ge_cortical_89pcw, aes(x = group,  y = HNRNPD, group = group, colour = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #axis.text.y = element_text(face="bold", color="black")
  ) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle= paste0(toString(gene), ' expression in ganglionic eminence and cortex at 8pcw and 9pcw'),
       x = 'structure',
       y = paste0(toString(gene),' expression (TPM)'))

h + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))

#CPSF7
gene <- gene_names[8] 

c <- ggplot(data = ge_cortical_89pcw, aes(x = group,  y = CPSF7, group = group, colour = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #axis.text.y = element_text(face="bold", color="black")
  ) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle= paste0(toString(gene), ' expression in ganglionic eminence and cortex at 8pcw and 9pcw'),
       x = 'structure',
       y = paste0(toString(gene),' expression (TPM)'))

c + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))


####################################################################################

#count how many samples are from each of the different brain regions in this dataset:
table(combined$structure_name)

#repeat the analysis for above genes, this time including more available (all cortex) samples at 8pcw and 9 pcw 
# to see if there is any difference

#keep only samples for 8pcw and 9pcw:
pcw_8_9 <- combined[combined$age %in% c("8 pcw", "9 pcw"),]

#keep only cortex and ganglionic eminence samples:
ge_cortex_all <- dplyr::filter(pcw_8_9, grepl('ganglionic eminence|cortex', structure_name))

#make a new column called 'region' labeling whether it is cortex or ganglionic eminences:
ge_cortex_all$region <- ifelse(grepl("GE", ge_cortex_all$structure_acronym), "ganglionic eminence", "cortex")


plotting(dataframe = ge_cortex_all, gene = "RUVBL2")
plotting(dataframe = ge_cortex_all, gene = "SOX2")
plotting(dataframe = ge_cortex_all, gene = "DCX")
plotting(dataframe = ge_cortex_all, gene = "CPSF7")
plotting(dataframe = ge_cortex_all, gene = "GABRA1")
plotting(dataframe = ge_cortex_all, gene = "PAX6")
plotting(dataframe = ge_cortex_all, gene = "MAP2")
plotting(dataframe = ge_cortex_all, gene = "SRSF9")
plotting(dataframe = ge_cortex_all, gene = "TCERG1")
plotting(dataframe = ge_cortex_all, gene = "HNRNPD")


#########################################################


#examine neuronal genetic markers in the BrainSpan RNA Seq dataset:
#NeuN has many aliases (different gene symbols) https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=146713
#CTIP2 also has many other gene symbols:  https://www.ncbi.nlm.nih.gov/gene/64919

markers <- TPM[TPM$gene_symbol %in% c("DCX", "GABRA1","GABA1", "CTIP2", "NeuN", "NEUN", "FOX3", "FOX-3", "HRNBP3","RBFOX3", 
                                      "ATL1", "RIT", "BCL11B","CTIP2", "IMD49", "CTIP-2","IDDFSTA", "SMARCM2",
                                      "ZNF856B", "ATL1-beta", "ATL1-alpha", "ATL1-delta", "ATL1-gamma", "hRIT1-alpha", "RIT1"),]
row.names(markers)<- markers$gene_symbol
markers <- markers[,2:525]

markers_t <- t(markers)
#merge markers expression dataset with the columns meta data
markers_meta <- cbind(markers_t, columns_meta)

View(markers_meta)

#filter for the samples which come from the cortex:
cortex <- dplyr::filter(markers_meta, grepl('cortex', structure_name))

View(cortex)

#plot the neuronal genetic markers expression data for cortical regions over all available time points

cortex$age <- as.factor(cortex$age)

#specify the order of ages:
cortex$age <- ordered(cortex$age, levels = c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw",
                                             "25 pcw", "26 pcw", "35 pcw", "37 pcw", "4 mos",  "10 mos", "1 yrs", "2 yrs","3 yrs",
                                             "4 yrs", "8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", "19 yrs", "21 yrs", "23 yrs", 
                                             "30 yrs", "36 yrs", "37 yrs","40 yrs"))
#GABRA1 in cortex all time points
ggplot(data = cortex, aes(x = structure_name,  y = GABRA1, colour = structure_name)) +
    geom_boxplot() +
    geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
    facet_wrap(~ age) +  #separate the boxplots by age
    theme_bw() + #set the background to white
    theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title = 'BrainSpan Transcriptome Dataset',
         subtitle=paste0('GABRA1 expression in cortex at 8pcw to 40 years'),
         x = 'structure',
         y = paste0('GABRA1',' expression (TPM)'))
  

plotting_cortex(dataframe = cortex, "DCX")
plotting_cortex(dataframe = cortex, "BCL11B")
plotting_cortex(dataframe = cortex, "RBFOX3")
colnames(cortex)

#############plot RUVBL2 expression over all prenatal time points (8pcw to 24 pcw) for all brain structures

#filter for prenatal time points (8pcw to 24 pcw):
prenatal <- combined[1:215,]

prenatal$age <- as.factor(prenatal$age)

#specify the order of ages:
prenatal$age <- ordered(prenatal$age, levels = c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw",
                                             "25 pcw"))

gene_name = "RUVBL2"

#NOTE: make sure the gene expression column is numeric, otherwise the y-axis labels will just be a blur together!!!
prenatal$RUVBL2<- as.numeric(prenatal$RUVBL2)

ggplot(data = prenatal, aes(x = structure_name,  y = RUVBL2, colour = structure_name)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle=paste0(gene_name, ' expression in all brain structures from 8pcw to 24 pcw'),
       x = '',
       y = paste0(gene_name,' expression (TPM)'))

# p +  theme(axis.text.y = element_text(face="bold", color="black", 
#                                       size=8),
#            axis.text.y = element_blank())
p +  theme(axis.text.y = element_blank())



##################### plot RUVBL2 expression in all available time points from 8pcw to 40 years for ganglionic eminences
# and all cortex regions

#convert age to factor
combined$age <- as.factor(combined$age)

#keep only cortex and ganglionic eminence samples:
ge_cortex_RUVBL2 <- dplyr::filter(combined, grepl('ganglionic eminence|cortex', structure_name))

#specify the order of ages:
ge_cortex_RUVBL2$age <- ordered(ge_cortex_RUVBL2$age, levels = c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw",
                                             "25 pcw", "26 pcw", "35 pcw", "37 pcw", "4 mos",  "10 mos", "1 yrs", "2 yrs","3 yrs",
                                             "4 yrs", "8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", "19 yrs", "21 yrs", "23 yrs", 
                                             "30 yrs", "36 yrs", "37 yrs","40 yrs"))

ggplot(data = ge_cortex_RUVBL2, aes(x = structure_name,  y = RUVBL2, colour = structure_name)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ age) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = 'BrainSpan Transcriptome Dataset',
       subtitle=paste0('RUVBL in ganglionic eminences and cortex over time'),
       x = 'structure',
       y = paste0('RUVBL2 expression (TPM)'))

################

# Save an object to a file
saveRDS(TPM, file = "TPM.rds")



####  Oct 20 2022  FABP7 in BrainSpan dataset over all ages; RNA expression in TPM 

#combined is the dataframe containing FABP7 expression, merged with columns_meta, then filtered for cortex and ganglionic eminences

#specify the order of ages:
FABP7$age <- ordered(FABP7$age, levels = c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw",
                                             "25 pcw", "26 pcw", "35 pcw", "37 pcw", "4 mos",  "10 mos", "1 yrs", "2 yrs","3 yrs",
                                             "4 yrs", "8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", "19 yrs", "21 yrs", "23 yrs", 
                                             "30 yrs", "36 yrs", "37 yrs","40 yrs"))




plotting_ge_cortex(dataframe = FABP7, "FABP7")

#make the same plot for FABP3, FABP5 and VRK1
#extract the rows for these genes in the TPM expression df
FABP <- TPM[TPM$gene_symbol %in% c("FABP3", "FABP5","VRK1"),]

#keep only expression rows
FABP<- FABP[,2:526]

FABP<-t(FABP)                            

#set the gene names as column names
colnames(FABP)<-FABP[1,]

#remove redundant gene_symbol row
FABP<-FABP[2:525,]

##remove the "X" from row names in data dataframe:
row.names(FABP) <- row.names(FABP) %>% str_replace("X","")

#combine RNA TPM expression with meta data
combined <- cbind(FABP, columns_meta)

#keep only cortex and ganglionic eminence samples:
combined <- dplyr::filter(combined, grepl('ganglionic eminence|cortex', structure_name))

#set age as factor
combined$age <- as.factor(combined$age)

#specify the order of ages:
combined$age <- ordered(combined$age, levels = c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw",
                                           "25 pcw", "26 pcw", "35 pcw", "37 pcw", "4 mos",  "10 mos", "1 yrs", "2 yrs","3 yrs",
                                           "4 yrs", "8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", "19 yrs", "21 yrs", "23 yrs", 
                                           "30 yrs", "36 yrs", "37 yrs","40 yrs"))

#convert the gene expression columns to numeric first to make sure the y-labels aren't jumbled together
combined$FABP3 <- as.numeric(combined$FABP3)
combined$FABP5 <- as.numeric(combined$FABP5)
combined$VRK1 <- as.numeric(combined$VRK1)

plotting_ge_cortex(dataframe = combined, "FABP3")
plotting_ge_cortex(dataframe = combined, "FABP5")
plotting_ge_cortex(dataframe = combined, "VRK1")

plotting_ge_cortex(dataframe = FABP7_8_9_pcw, "FABP7")

############### February 2 2023

#comparing FABP7 in cortex versus ganglionic eminences at 17 pcw

#subset for 17 pcw only
FABP7_17pcw <- FABP7[FABP7$age == "17 pcw", ]

FABP7_16pcw <- FABP7[FABP7$age == "16 pcw", ]

both<- rbind(FABP7_16pcw, FABP7_17pcw)

#keep only FABP7 expression levels, age and structure name columns
both <- both[,c(1, 5, 9)]



setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Allen Brain Atlas/BrainSpan/data/")

#extract 16pcw and 17pcw FABP7 data into csv file 
write.csv(both, "16pc_17pcw_FABP7.csv")

#manually add in NA where there are no measurements
FABP7_pcw_17pcw_NAs<-read.csv("16pc_17pcw_FABP7.csv")

FABP7_pcw_17pcw_NAs<-FABP7_pcw_17pcw_NAs[,2:4]

table(FABP7_pcw_17pcw_NAs$structure_name)

#rearrange order of columns

FABP7 <- cbind(FABP7_pcw_17pcw_NAs$structure_name,FABP7_pcw_17pcw_NAs$age)

FABP7<- cbind(FABP7, FABP7_pcw_17pcw_NAs$FABP7)

colnames(FABP7)<-c("structure_name", "age", "FABP7")

FABP7df<- as.data.frame(FABP7)

FABP7df$FABP7 <- as.numeric(FABP7df$FABP7)
FABP7df$age<-as.character(FABP7df$age)

#assign a numerical ID to the structure names to see if it works with the imputations 
FABP7df$structure_ID <- nrow(FABP7df)
FABP7df$structure_ID<-rep(1:13,each=6)

FABP7_structure_key<- cbind(FABP7df[,1],FABP7df[,4])

FABP7df <- cbind(FABP7df[,4], FABP7df[,2:3])
colnames(FABP7df)[1]<- "structure_ID"

#do a repeated measure test on FABP7 expression between 16pcw and 17 pcw 
#impute NA values
ini <- mice(data = FABP7df, maxit = 0)


#define the 'subject' identifier (code as '-2' in predictor matrix)
pred <- ini$pred
pred["FABP7", "structure_ID"] <- -2

#run MI
#"2l.pan is level-1 normal homoscedastic, lmer
imp <- mice(data = FABP7df, pred = pred, method = "2l.pan", m = 100)

summary(imp)

#create a list of completed data sets
implist <- mids2mitml.list(imp)

#fit the ANOVA model using lmer(). The model contains an additional term (1|electrode_number), which
#specifies a random effect for each subject. This effect captures unsystematic differences between subjects,
#thus accounting for the nested structure of the repeated-measures data.

#fit the ANOVA model
fit3 <- with(implist, lmer(FABP7 ~ 1 + age + (1|structure_ID)))
testEstimates(fit3, var.comp = TRUE)

# Call:
#   
#   testEstimates(model = fit3, var.comp = TRUE)
# 
# Final parameter estimates and inferences obtained from 100 imputed data sets.
# 
#               Estimate Std.Error   t.value        df   P(>|t|)       RIV       FMI 
# (Intercept)   830.097    88.980     9.329  2441.846     0.000     0.252     0.202 
# age17 pcw      89.582   125.300     0.715   648.746     0.475     0.641     0.393 
# 
# Estimate 
# Intercept~~Intercept|structure_ID 2.002e+04 
# Residual~~Residual                1.866e+05 
# ICC|structure_ID                  9.895e-02 
# 
# Unadjusted hypothesis test as appropriate in larger samples.
# 
# Warning message:
#   In .checkDeprecated(extra.pars, arg.list = dots, name = "var.comp") :
#   The 'var.comp' argument is deprecated. Please use 'extra.pars' instead.

#pool the parameter estimates
fit3.reduced <- with(implist, lmer(FABP7 ~ 1 + (1|structure_ID)))

testModels(fit3, fit3.reduced, method = "D1")

# Call:
#   
#   testModels(model = fit3, null.model = fit3.reduced, method = "D1")
# 
# Model comparison calculated from 100 imputed data sets.
# Combination method: D1
# 
# F.value     df1     df2   P(>F)     RIV 
# 0.511       1   611.301   0.475   0.641 
# 
# Unadjusted hypothesis test as appropriate in larger samples.
# Models fitted with REML were used as is.

### Feb 7, 2023
#Compare everything to an early time point when FABP7 expression is low

#like 8 pcw vs X....
#since 16 and 17 look similar so not suprising -p value not super low


#8 pcw vs. 9 pcw
fabp7_8_9pcw<-get_tp(combined_df = combined, timepoint1 = "8 pcw", timepoint2 = "9 pcw")

#sort by structure_name
write.csv(fabp7_8_9pcw, "fabp7_8_9pcw.csv")

#put in NA manually  for the structures which are missing FABP7 values at 1 of the 2 timepoints
fabp7_8_9pcw<-read.csv("fabp7_8_9pcw.csv")

#keep just FABP7, age and structure name
fabp7_8_9pcw <-fabp7_8_9pcw[,c(1,2,4)]

#assign a numerical ID to the structure names to work with the imputations 
fabp7_8_9pcw$structure_ID <- nrow(fabp7_8_9pcw)

#count how many unique structures there are in this comparison
length(unique(fabp7_8_9pcw$structure_name))
fabp7_8_9pcw$structure_ID<-rep(1:17,each=2)

#keep structure ID key and structure name in a dataframe for reference
fabp7_8_9pcw_structure_key <-fabp7_8_9pcw[,c(3,4)]
FABP7_8_9pcw <- cbind(fabp7_8_9pcw[,4], fabp7_8_9pcw[,c(1,2)])
colnames(FABP7_8_9pcw)[1]<- "structure_ID"

#impute NA values
ini <- mice(data = FABP7_8_9pcw, maxit = 0)

#define the 'subject' identifier (code as '-2' in predictor matrix)
pred <- ini$pred
pred["FABP7", "structure_ID"] <- -2

#run MI
#"2l.pan is level-1 normal homoscedastic, lmer
imp <- mice(data = FABP7_8_9pcw, pred = pred, method = "2l.pan", m = 100)

summary(imp)

#create a list of completed data sets
implist <- mids2mitml.list(imp)

#fit the ANOVA model using lmer(). The model contains an additional term (1|electrode_number), which
#specifies a random effect for each subject. This effect captures unsystematic differences between subjects,
#thus accounting for the nested structure of the repeated-measures data.

#fit the ANOVA model
fit3 <- with(implist, lmer(FABP7 ~ 1 + age + (1|structure_ID)))
testEstimates(fit3, var.comp = TRUE)

# Call:
#   
#   testEstimates(model = fit3, var.comp = TRUE)
# 
# Final parameter estimates and inferences obtained from 100 imputed data sets.
# 
# Estimate Std.Error   t.value        df   P(>|t|)       RIV       FMI 
# (Intercept)   365.929    54.827     6.674 32317.443     0.000     0.059     0.055 
# age9 pcw      -24.739    63.000    -0.393  2039.628     0.695     0.283     0.221 
# 
# Estimate 
# Intercept~~Intercept|structure_ID 21969.123 
# Residual~~Residual                26303.729 
# ICC|structure_ID                      0.458 
# 
# Unadjusted hypothesis test as appropriate in larger samples.

#pool the parameter estimates
fit3.reduced <- with(implist, lmer(FABP7 ~ 1 + (1|structure_ID)))

testModels(fit3, fit3.reduced, method = "D1")
# 
# Call:
#   
#   testModels(model = fit3, null.model = fit3.reduced, method = "D1")
# 
# Model comparison calculated from 100 imputed data sets.
# Combination method: D1
# 
# F.value      df1      df2    P(>F)      RIV 
# 0.154        1 1900.047    0.695    0.283 
# 
# Unadjusted hypothesis test as appropriate in larger samples.
# Models fitted with REML were used as is.


################################### 8 pcw vs. 12 pcw
#test for differences in FABP7 expression between ganglionic eminences & cortex regions at 8 pcw versus
#FABP7 expression in ganglionic eminences and cortex regions at 12 pcw
fabp7_8_12pcw<-get_tp(combined_df = FABP7, timepoint1 = "8 pcw", timepoint2 = "12 pcw")

table(fabp7_8_12pcw$structure_name)
#keep only relevant columns
fabp7_8_12pcw_filtered<-fabp7_8_12pcw[,c(1, 5, 9)]

#sort by Structure name, then fill in NA values for missing values
write.csv(fabp7_8_12pcw_filtered, "fabp7_8_12pcw.csv")

fabp7_8_12pcw <- read.csv("fabp7_8_12pcw.csv")

#keep just FABP7, age and structure name

#assign a numerical ID to the structure names to work with the imputations 
fabp7_8_12pcw$structure_ID <- nrow(fabp7_8_12pcw)

#count how many unique structures there are in this comparison
length(unique(fabp7_8_12pcw$structure_name))
fabp7_8_12pcw$structure_ID<-rep(1:18,each=6)

#keep structure ID key and structure name in a dataframe for reference
fabp7_8_12pcw_structure_key <-fabp7_8_12pcw[,c(3,4)]
FABP7_8_12pcw <- cbind(fabp7_8_12pcw[,4], fabp7_8_12pcw[,c(1,2)])
colnames(FABP7_8_12pcw)[1]<- "structure_ID"

fabp7_8_12pcw_imputed<- impute_process(test_df = FABP7_8_12pcw)
fabp7_8_12pcw_imputed
# Call:
#   
#   testEstimates(model = fit3, var.comp = TRUE)
# 
# Final parameter estimates and inferences obtained from 100 imputed data sets.
# 
#               Estimate Std.Error   t.value        df   P(>|t|)       RIV       FMI 
# (Intercept)   225.802    36.873     6.124   704.485     0.000     0.600     0.377 
# age8 pcw       80.473    50.385     1.597   755.176     0.111     0.568     0.364 
# 
# Estimate 
# Intercept~~Intercept|structure_ID   723.097 
# Residual~~Residual                43726.425 
# ICC|structure_ID                      0.016 
# 
# Unadjusted hypothesis test as appropriate in larger samples.

fit3.reduced_fabp7_8_12pcw <-fit_anova(implist)
fit3.reduced_fabp7_8_12pcw 

# Call:
#   
#   testModels(model = fit3, null.model = fit3.reduced, method = "D1")
# 
# Model comparison calculated from 100 imputed data sets.
# Combination method: D1
# 
# F.value     df1     df2   P(>F)     RIV 
#   2.551       1 710.106   0.111   0.568 
# 
# Unadjusted hypothesis test as appropriate in larger samples.
# Models fitted with REML were used as is.



################################### 8 pcw vs. 16 pcw
fabp7_8_16pcw<-get_tp(combined_df = FABP7, timepoint1 = "8 pcw", timepoint2 = "16 pcw")

write.csv(fabp7_8_16pcw, "fabp7_8_16pcw.csv")

fabp7_8_16pcw_NA<- read.csv("fabp7_8_16pcw.csv")
fabp7_8_16pcw_NA$structure_ID<-rep(1:18,each=6)

#keep structure ID key and structure name in a dataframe for reference
fabp7_8_16pcw_structure_key <-fabp7_8_16pcw_NA[,c(3,4)]
FABP7_8_16pcw <- cbind(fabp7_8_16pcw_NA[,4], fabp7_8_16pcw_NA[,c(1,2)])
colnames(FABP7_8_16pcw)[1]<- "structure_ID"

fabp7_8_16pcw_imputed<- impute_process(test_df = FABP7_8_16pcw)
fabp7_8_16pcw_imputed

# Call:
#   
#   testEstimates(model = fit3, var.comp = TRUE)
# 
# Final parameter estimates and inferences obtained from 100 imputed data sets.
# 
#               Estimate Std.Error   t.value        df   P(>|t|)       RIV       FMI 
# (Intercept)   761.400    75.560    10.077   550.454     0.000     0.736     0.426 
# age8 pcw     -116.966    96.503    -1.212   724.644     0.226     0.586     0.371 
# 
# Estimate 
# Intercept~~Intercept|structure_ID 6.348e+03 
# Residual~~Residual                1.585e+05 
# ICC|structure_ID                  3.911e-02 
# 
# Unadjusted hypothesis test as appropriate in larger samples.

fit3.reduced_fabp7_8_16pcw <-fit_anova(implist)
fit3.reduced_fabp7_8_16pcw 

# Call:
#   
#   testModels(model = fit3, null.model = fit3.reduced, method = "D1")
# 
# Model comparison calculated from 100 imputed data sets.
# Combination method: D1
# 
# F.value     df1     df2   P(>F)     RIV 
# 1.469       1 681.767   0.226   0.586 
# 
# Unadjusted hypothesis test as appropriate in larger samples.
# Models fitted with REML were used as is.


######################### 8 pcw vs 17 pcw...
fabp7_8_17pcw<-get_tp(combined_df = FABP7, timepoint1 = "8 pcw", timepoint2 = "17 pcw")

write.csv(fabp7_8_17pcw, "fabp7_8_17pcw.csv")
#add NAs where there's missing values

fabp7_8_17pcw_NA <-read.csv("fabp7_8_17pcw.csv")
fabp7_8_17pcw_NA$structure_ID<-rep(1:16,each=2)

#keep structure ID key and structure name in a dataframe for reference
fabp7_8_17pcw_structure_key <-fabp7_8_17pcw_NA[,c(3,4)]
FABP7_8_17pcw <- cbind(fabp7_8_17pcw_NA[,4], fabp7_8_17pcw_NA[,c(1,2)])
colnames(FABP7_8_17pcw)[1]<- "structure_ID"

fabp7_8_17pcw_imputed<- impute_process(test_df = FABP7_8_17pcw)
fabp7_8_17pcw_imputed

# Call:
#   
#   testEstimates(model = fit3, var.comp = TRUE)
# 
# Final parameter estimates and inferences obtained from 100 imputed data sets.
# 
#               Estimate Std.Error   t.value        df   P(>|t|)       RIV       FMI 
# (Intercept)   936.295   154.138     6.074   629.211     0.000     0.657     0.399 
# age8 pcw     -424.578   189.879    -2.236   762.977     0.026     0.563     0.362 
# 
# Estimate 
# Intercept~~Intercept|structure_ID 4.482e+04 
# Residual~~Residual                1.845e+05 
# ICC|structure_ID                  1.938e-01 
# 
# Unadjusted hypothesis test as appropriate in larger samples.

fit3.reduced_fabp7_8_17pcw <-fit_anova(implist)
fit3.reduced_fabp7_8_17pcw 

# Call:
#   
#   testModels(model = fit3, null.model = fit3.reduced, method = "D1")
# 
# Model comparison calculated from 100 imputed data sets.
# Combination method: D1
# 
# F.value     df1     df2   P(>F)     RIV 
#  5.000       1 717.346   0.026   0.563 
# 
# Unadjusted hypothesis test as appropriate in larger samples.
# Models fitted with REML were used as is.