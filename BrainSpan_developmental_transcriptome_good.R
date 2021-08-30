#Allen Human Brain Atlas Developmental Transcriptome RNA Seq data Analysis 

#Data downloaded from https://www.brainspan.org/static/download.html RNA-Seq Gencode v10 summarized to genes

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

#set working directory to data directory 
# data <- read.table("expression_matrix.txt", header=T, row.names=1, sep ="\t")
data <- read.csv(file = 'expression_matrix.csv', header=F)  #expression values in RPKM
columns_meta <- read.csv(file = 'columns_metadata.csv', header = T)
rows_meta <- read.csv(file = 'rows_metadata.csv', header = T)

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

###########################read in RNA expression data in TPM
tpm <- readRDS("TPM.rds")
columns_meta <- read.csv(file = 'columns_metadata.csv', header = T)

#filter for genes of interest:
genes <- tpm[tpm$gene_symbol %in% c("RUVBL2", "SOX2","GABRA1", "PAX6", "MAP2", "DCX", "SRSF9", "TCERG1", "HNRNPD", "CPSF7",
                                      "GABRA1","GABA1", "CTIP2", "NeuN", "NEUN", "FOX3", "FOX-3", "HRNBP3","RBFOX3", 
                                       "ATL1", "RIT", "BCL11B","CTIP2", "IMD49", "CTIP-2","IDDFSTA", "SMARCM2",
                                       "ZNF856B", "ATL1-beta", "ATL1-alpha", "ATL1-delta", "ATL1-gamma", "hRIT1-alpha", "RIT1"),]

#merge RNA expression data (TPM) with columns meta data:
row.names(genes)<- genes$gene_symbol
genes$gene_symbol <- NULL
genes_t <- t(genes)

combined <- cbind(genes_t, columns_meta)

combined$age <- as.factor(combined$age)

#specify the order of ages:
combined$age <- ordered(combined$age, levels = c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw",
                                             "25 pcw", "26 pcw", "35 pcw", "37 pcw", "4 mos",  "10 mos", "1 yrs", "2 yrs","3 yrs",
                                             "4 yrs", "8 yrs", "11 yrs", "13 yrs", "15 yrs", "18 yrs", "19 yrs", "21 yrs", "23 yrs", 
                                             "30 yrs", "36 yrs", "37 yrs","40 yrs"))

#keep only cortex and ganglionic eminence samples:
ge_cortex <- dplyr::filter(combined, grepl('ganglionic eminence|cortex', structure_name))

#filter for prenatal time points (8pcw to 24 pcw):
prenatal <- ge_cortex[ge_cortex$age %in% c("8 pcw", "9 pcw",  "12 pcw", "13 pcw", "16 pcw", "17 pcw", "19 pcw", "21 pcw", "24 pcw"),]

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
       subtitle=paste0('RUVBL in ganglionic eminences and cortex during prenatal development'),
       x = 'structure',
       y = paste0('RUVBL2 expression (TPM)'))
