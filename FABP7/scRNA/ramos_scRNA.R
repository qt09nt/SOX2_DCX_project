#Analysis of FABP7 expression in  "Late prenatal human development single-nucleus 
#transcriptomics" data from Ramos et al (2022) paper

#single cell data downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511

#January 19 2023

library(Seurat)

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/late prenatal human neurodevelopment resolved by single nucleus transcriptomics/processed/")

#test <- Matrix::readMM('matrix.mtx')

#load in processed data in TPM format
annotation <-read.table(file = 'Human.GRCh38.p13.annot.tsv', sep = '\t', header = TRUE, fill=TRUE)
data <- read.table(file = 'GSE217511_norm_counts_tpm_ncbi.tsv', sep = '\t', header = TRUE, fill=TRUE)

#try re-reading the annotation and data files which have been saved as tab delimited text files
#and where only the row containing FABP7 expression values are kept
annotation <-read.table(file = 'Human.GRCh38.p13.annot_FABP7_only.txt', sep = '\t', header = TRUE)
data <- read.table(file = 'GSE217511_norm_counts_tpm_ncbi_FABP7_only.txt', sep = '\t', header = TRUE)
meta <- read.csv(file = "meta.csv")

meta$X <- NULL

#specify the order of the ages
meta$developmental.stage <- ordered(meta$developmental.stage, levels = c("17 weeks gestation",   "22 weeks gestation",   "23.3 weeks gestation", "23.7 weeks gestation",
                                                                         "24.7 weeks gestation", "26 weeks gestation",   "32.8 weeks gestation", "33.2 weeks gestation",
                                                                         "38.3 weeks gestation", "28 years",             "45 years",             "53 years"        ))

data_t <- t(data)
colnames(data_t)[1]<-"FABP7"
row.names(meta) <- meta$sample

#combine the expression data df and meta df
fabp7<- cbind(data_t, meta)

ggplot(data = fabp7, aes(x = tissue,  y = FABP7, group = tissue, colour = tissue)) +
  geom_boxplot() +
  geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
  facet_wrap(~ developmental.stage) +  #separate the boxplots by age
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  ) +
  labs(title = 'FABP7 Expression in Ramos et al (2022) Transcriptome Dataset',
       subtitle='',
       x = 'tissue',
       y = 'FABP7 expression - RNA Seq (TPM)')


#Read10X: Load in data from 10X
#In Seurat: Tools for Single Cell Genomics 
#https://rdrr.io/cran/Seurat/man/Read10X.html






















#remove duplicated Gene IDs
annotation_no_dup = annotation[!duplicated(annotation$GeneID),]
row.names(annotation_no_dup)<-annotation_no_dup$GeneID

data_no_dup = data[!duplicated(data$GeneID),]
row.names(data_no_dup)<-data_no_dup$GeneID


#merge the annotation file and data file
combined <- merge(annotation, data,
                        by = 'row.names', all = TRUE)
