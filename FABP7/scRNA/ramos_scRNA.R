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



#January 24 2023


####try to load in the raw files

# For output from CellRanger >= 3.0 with multiple data types
data_dir <- 'path/to/data/directory'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)

data_dir <- 'C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/late prenatal human neurodevelopment resolved by single nucleus transcriptomics/raw/GSM6720852_9C'
list.files(data_dir)
data <- Read10X(data.dir = data_dir)
# Error in if (any(els$j < 1 | els$j > nc)) stop("readMM(): column values 'j' are not in 1:nc",  : 
# missing value where TRUE/FALSE needed
# In addition: Warning message:
# In scan(file, nmax = nz, quiet = TRUE, what = list(i = integer(),  :
# number of items read is not a multiple of the number of columns
seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)

#try another sample
data_dir <- 'C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/late prenatal human neurodevelopment resolved by single nucleus transcriptomics/raw/GSM6720853_9G'
list.files(data_dir)
data <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)



############ Seurat tutorial 
#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(dplyr)
library(Seurat)
library(patchwork)

#Load data from Ramos et al 36 samples raw files
setwd("D:/Technical Analyst I UHN May 4 2021/organoid group/Sofia/raw/GSM6720852")
GSM6720852.data <- Read10X(data.dir = "D:/raw/GSM6720852/")
GSM6720852 <- CreateSeuratObject(counts = GSM6720852.data, project = 'ramos', min.cells =3, min.features = 200)
GSM6720852
# An object of class Seurat 
# 25066 features across 58897 samples within 1 assay 
# Active assay: RNA (25066 features, 0 variable features)

# Standard pre-processing workflow
# 
# The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. These represent the selection and filtration of cells based on QC metrics, data normalization and scaling, and the detection of highly variable features.
# QC and selecting cells for further analysis
# 
# Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics commonly used by the community include
# 
# The number of unique genes detected in each cell.
# Low-quality cells or empty droplets will often have very few genes
# Cell doublets or multiplets may exhibit an aberrantly high gene count
# Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
# The percentage of reads that map to the mitochondrial genome
# Low-quality / dying cells often exhibit extensive mitochondrial contamination
# We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
# We use the set of all genes starting with MT- as a set of mitochondrial genes


#the [[ ]] operator can add columns to object metadata. This is a great place to stash QC stats
GSM6720852[["percent.mt"]]<- PercentageFeatureSet(GSM6720852, pattern = "^MT-")

#show QC metrics for the first 5 cells
head(GSM6720852@meta.data, 5)

#                    orig.ident nCount_RNA nFeature_RNA percent.mt
# AAACCCAAGAATTGCA-1      ramos        472          409  2.5423729
# AAACCCAAGACAGCGT-1      ramos        546          447  3.4798535
# AAACCCAAGACCCTTA-1      ramos        509          428  0.0000000
# AAACCCAAGAGGCGGA-1      ramos        557          454  0.5385996
# AAACCCAAGATGACCG-1      ramos        554          449  1.6245487

###########Visualize QC metrics and use these to filter cells

#filter cells that have unique feature counts over 2500 or less than 200
#filter cells that have >5% mitochondrial counts

#Visualize QC metrics as a violin plot
VlnPlot(GSM6720852, features = c("nFeature_RNA","nCount_RNA", "percent.mt"), ncol = 3)

