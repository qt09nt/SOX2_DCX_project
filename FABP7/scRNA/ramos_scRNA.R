#Analysis of FABP7 expression in  "Late prenatal human development single-nucleus 
#transcriptomics" data from Ramos et al (2022) paper. This is the analysis done on the pre-processed data which is available in the Supplementary section of 
#Ramos et al (2022) paper.

#For these plots I worked with this file Combined_log-normalized_AveExp_snRNAseq_Human.xlsx found at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6720852: 
#The file contains processed, log-normalized average gene expression data calculated for each of the following integrated objects: all prenatal germinal 
#matrix samples, subclustered glia and progenitors from all prenatal germinal matrix samples, all prenatal cortical plate samples, subclustered glia and 
#progenitors from 24-41gw samples, all adult SVZ+caudate samples, all adult neocortex samples. All clusters are also annotated by cell type / state. 
#UMAPS of each object are included along with cell type abbreviations.

#single cell data downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217511

#January 19 2023

library(Seurat)
library(umap)
library(dplyr)
library(cowplot)
library(ggplot2)
library(M3C)  
library(stringr)
library(patchwork)

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

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(GSM6720852, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(GSM6720852, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

GSM6720852 <- subset(GSM6720852, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)





# January 26, 2023
################################ working with the supplemental Combined_log-normalized_AveExp_snRNAseq_Human.xlsx file 
# downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6720852

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/late prenatal human neurodevelopment resolved by single nucleus transcriptomics/supplemental")

abbreviations <- read.csv("abbreviations.csv")

#load data for Processed, log-normalized average gene expression data calculated for
#each of the following integrated objects: all prenatal germinal matrix samples
prenatal_gm_all <-read.csv("prenatal_gm_all.csv")
prenatal_cp_all <-read.csv("prenatal_cp_all.csv")
adult_svz_caudate <- read.csv("adult_svz_caudate.csv")
adult_neocortex <-read.csv("adult_neocortex.csv")
prenatal_gm_subclustered<-read.csv("prenatal_gm_subclustered.csv")
prenatal_cp_subclust <- read.csv("prenatal_cp_subclustered.csv")

colnames(prenatal_gm_all)[1]<-"gene"
colnames(prenatal_cp_all)[1]<-"gene"
colnames(adult_svz_caudate)[1]<- "gene"
colnames(adult_neocortex)[1]<-"gene"
colnames(prenatal_gm_subclustered)[1]<-"gene"
colnames(prenatal_cp_subclust)[1]<-"gene"

colnames(abbreviations)[1]<-"abbreviation"

#split the abbreviations column into "abbreviation" and "name"
abbreviations[c("abbrev","name")]<- str_split_fixed(abbreviations$abbreviation, " = ", 2)

#https://cran.r-project.org/web/packages/umap/vignettes/umap.html


#remove duplicated gene IDs 
prenatal_gm_all = prenatal_gm_all[!duplicated(prenatal_gm_all$gene),]
prenatal_cp_all = prenatal_cp_all[!duplicated(prenatal_cp_all$gene),]
adult_svz_caudate = adult_svz_caudate[!duplicated(adult_svz_caudate$gene),]
adult_neocortex = adult_neocortex[!duplicated(adult_neocortex$gene),]
prenatal_gm_subclustered = prenatal_gm_subclustered[!duplicated(prenatal_gm_subclustered$gene),]
prenatal_cp_subclust = prenatal_cp_subclust[!duplicated(prenatal_cp_subclust$gene),]

row.names(prenatal_gm_all)<-prenatal_gm_all$gene
row.names(prenatal_cp_all)<-prenatal_cp_all$gene
row.names(adult_svz_caudate)<-adult_svz_caudate$gene
row.names(adult_neocortex)<-adult_neocortex$gene
row.names(prenatal_gm_subclustered)<-prenatal_gm_subclustered$gene
row.names(prenatal_cp_subclust)<-prenatal_cp_subclust$gene

#remove redundant gene column
prenatal_gm_all$gene <- NULL
prenatal_cp_all$gene <- NULL
adult_svz_caudate$gene <- NULL
adult_neocortex$gene<- NULL
prenatal_gm_subclustered$gene<- NULL
prenatal_cp_subclust$gene<-NULL

#get cell types of the normalized samples
cell_type<- colnames(prenatal_gm_all)
cell_type_prenatal_gm_all <- as.data.frame(cell_type)

cell_type<-colnames(prenatal_cp_all)
cell_type_prenatal_cp_all<- as.data.frame(cell_type)

cell_type<-colnames(adult_svz_caudate)
cell_type_adult_svz_caudate<- as.data.frame(cell_type)

cell_type<-colnames(adult_neocortex)
cell_type_adult_neocortex<- as.data.frame(cell_type)

cell_type<-colnames(prenatal_gm_subclustered)
cell_type_prenatal_gm_subclustered<- as.data.frame(cell_type)

cell_type<-colnames(prenatal_cp_subclust)
cell_type_prenatal_cp_subclust<- as.data.frame(cell_type)

#keep just the cell type part of the name
cell_type[,1]<-gsub(".*\\.","", cell_type[,1])

#label the cell type of the samples with the full name cell type labels
cell_type$cell_type_name <- "NA"


#function for processing the cell type dataframe
process_celltype <-function(cell_type_df){
  
  cell_type<-as.data.frame(cell_type_df)
  
  #keep just the cell type part of the name
  cell_type[,1]<-gsub(".*\\.","", cell_type[,1])
  
  #label the cell type of the samples with the full name cell type labels
  cell_type$cell_type_name <- "NA"
  
  
  #fill in cell type name based on condition
  cell_type$cell_type_name[cell_type$cell_type == "AC"] <- "Astrocyte"
  cell_type$cell_type_name[cell_type$cell_type == "BVC"] <- "Blood vessel cell (endothelium+pericytes)"
  cell_type$cell_type_name[cell_type$cell_type == "CP"] <- "Choroid plexus"
  cell_type$cell_type_name[cell_type$cell_type == "CPN"] <- "Cortical projection neuron"
  cell_type$cell_type_name[cell_type$cell_type == "gIPC"] <- "Glial intermediate progenitor cell"
  cell_type$cell_type_name[cell_type$cell_type == "IN"] <- "Interneuron"
  cell_type$cell_type_name[cell_type$cell_type == "MG"] <- "Microglila"
  cell_type$cell_type_name[cell_type$cell_type == "MSN"] <- "Medium spiny neuron (D1 or D2)"
  cell_type$cell_type_name[cell_type$cell_type == "nIPC"] <- "Neuronal intermediate progenitor cell"
  cell_type$cell_type_name[cell_type$cell_type == "OPC"] <- "Oligodendrocyte progenitor/precursor cell"
  cell_type$cell_type_name[cell_type$cell_type == "SPN"] <- "Subplate neuron"
  cell_type$cell_type_name[cell_type$cell_type == "TAC"] <- "transit amplifying cell / cycling progenitor"
  cell_type$cell_type_name[cell_type$cell_type == "UD"] <- "Undefined cell type or lineage"
  cell_type$cell_type_name[cell_type$cell_type == "EN"] <- "Excitatory neuron"
  cell_type$cell_type_name[cell_type$cell_type == "CRN"] <- "Cajal Retzius cell"
  cell_type$cell_type_name[cell_type$cell_type == "OL"] <-"Oligodendrocyte"
  cell_type$cell_type_name[cell_type$cell_type == "EPD"] <-"Ependymal cell"
  cell_type$cell_type_name[cell_type$cell_type == "p"] <-"Astrocyte, protoplasmic type"
  cell_type$cell_type_name[cell_type$cell_type == "pv"] <-"Interneuron.pv"
  cell_type$cell_type_name[cell_type$cell_type == "5ht"] <-"Interneuron.5ht"
  cell_type$cell_type_name[cell_type$cell_type == "f"] <-"Astrocyte, fibrous type"
  cell_type$cell_type_name[cell_type$cell_type == "som"] <-"Interneuron.som"
  cell_type$cell_type_name[cell_type$cell_type == "Tcell"] <- "T-lymphocyte"
  cell_type$cell_type_name[cell_type$cell_type =="preOL"] <-"Premyelinating / early myelinating BCAS1+ oligodendrocyte"
  cell_type$cell_type_name[cell_type$cell_type =="oRG"] <- "Outer radial glia"
  cell_type$cell_type_name[cell_type$cell_type == "mIPC"]<- "multipotent intermediate progenitor cell"
  cell_type$cell_type_name[cell_type$cell_type == "O"]<- "Glial intermediate progenitor cell, OPC bias"
  cell_type$cell_type_name[cell_type$cell_type == "A"]<-"Glial intermediate progenitor cell, Astrocyte bias"
  cell_type$cell_type_name[cell_type$cell_type == "tRG"]<-"Truncated radial glia"
  
  #write cell type to global env
  cell_type<<-cell_type
}

cell_type_adult_neocortex<-process_celltype(cell_type_adult_neocortex)
cell_type_adult_svz_caudate<-process_celltype(cell_type_adult_svz_caudate)
cell_type_prenatal_cp_all<-process_celltype(cell_type_prenatal_cp_all)
cell_type_prenatal_cp_subclust<-process_celltype(cell_type_prenatal_cp_subclust)
cell_type_prenatal_gm_all<-process_celltype(cell_type_prenatal_gm_all)
cell_type_prenatal_gm_subclustered<-process_celltype(cell_type_prenatal_gm_subclustered)


table(cell_type$cell_type)
# AC  BVC   CP  CPN gIPC   IN   MG  MSN nIPC  OPC  SPN  TAC   UD 
# 4    1    1   11    1    7    1    2    2    1    3    1    3 



cell_type$cell_type_name<- as.factor(cell_type$cell_type_name)

#keep complete cases only (ie. omit the row with the NA)
prenatal_gm_all_full <- prenatal_gm_all[complete.cases(prenatal_gm_all),]
prenatal_cp_all_full <- prenatal_cp_all[complete.cases(prenatal_cp_all),]
adult_svz_caudate_full <- adult_svz_caudate[complete.cases(adult_svz_caudate),]
adult_neocortex_full <- adult_neocortex[complete.cases(adult_neocortex),]
prenatal_gm_subclust_full <- prenatal_gm_subclustered[complete.cases(prenatal_gm_subclustered),]
prenatal_cp_subclust_full <- prenatal_cp_subclust[complete.cases(prenatal_cp_subclust),]

#since prenatal_cp_subclust only has 14 samples, n_neighbors should be changed to a number less than 14
custom.config<-umap.defaults
custom.config$n_neighbors <- 13
prenatal_cp_subclust_t<-t(prenatal_cp_subclust_full)
prenatal_cp_subclust.umap<-umap(prenatal_cp_subclust_t, config = custom.config)

#create list object with the prenatal_gm_all expression data and the cell type labels
prenatal_gm_all_combined <- list(expr = prenatal_gm_all, cellnames = cell_type$cell_type_name)

prenatal_cp_all_combined <- list(expr = prenatal_cp_all_full, cellnames = cell_type$cell_type_name)
adult_svz_caudate_combined <- list(expr = adult_svz_caudate_full, cellnames = cell_type$cell_type_name)

adult_neocortex_combined <- list(expr = adult_neocortex_full, cellnames = cell_type$cell_type_name)

prenatal_gm_subclust_combined <- list(expr = prenatal_gm_subclust_full, cellnames = cell_type$cell_type_name)

#prenatal_cp_subclust_combined <-list(expr = prenatal_cp_subclust_full, cellnames = cell_type$cell_type_name)

as.matrix(names(prenatal_gm_all)[colSums(is.na(prenatal_gm_all))!=0]) 
# [,1]             
# [1,] "RNA.24.L2.3.CPN"

#apparently there is a NA in this column

sum(is.na(prenatal_gm_all$RNA.24.L2.3.CPN))
#1



#create list object with the prenatal_gm_all expression data and the cell type labels
prenatal_gm_all_combined <- list(expr = prenatal_gm_all_full, cellnames = cell_type$cell_type_name)

library(umap)

#note if M3C not loaded, umap default library will just run UMAP without plotting
umap(pollen$data,labels=as.factor(pollen$celltypes),controlscale=TRUE,scale=3)

umap(prenatal_gm_all_combined$expr, labels=as.factor(prenatal_gm_all_combined$cellnames), controlscale=TRUE, scale=3)

prenatal_gm_all.umap <-umap(prenatal_gm_all_combined$expr, labels=as.factor(prenatal_gm_all_combined$cellnames), controlscale=TRUE, scale=3)

prenatal_cp_all.umap <- umap(prenatal_cp_all_combined$expr, labels = as.factor(prenatal_cp_all_combined$cellnames), controlscale = TRUE, scale =3)
adult_svz_caudate.umap <- umap(adult_svz_caudate_combined$expr, labels = as.factor(adult_svz_caudate_combined$cellnames), controlscale = TRUE, scale =3 )

adult_neocortex.umap <- umap(adult_neocortex_combined$expr, labels = as.factor(adult_neocortex_combined$cellnames), controlscale = TRUE, scale =3 )

prenatal_gm_subclust.umap<-umap(prenatal_gm_subclust_combined$expr, labels = as.factor(prenatal_gm_subclust_combined$cellnames), controlscale = TRUE, scale =3)


prenatal_cp_subclust.umap<-umap(prenatal_cp_subclust_combined$expr, config=custom.config, labels = as.factor(prenatal_cp_subclust_combined$cellnames), controlscale = TRUE, scale =3)

#prenatal_cp_subclust.umap<-umap(prenatal_cp_subclust_combined$expr, labels = as.factor(prenatal_cp_subclust_combined$cellnames), controlscale = TRUE, scale =3)


#get umap coordinates into a df
prenatal_cp_all_umap_coordinates<- prenatal_cp_all.umap[["data"]]
row.names(cell_type)<-row.names(prenatal_cp_all_umap_coordinates)

adult_svz_caudate_umap_coordinates <- adult_svz_caudate.umap[["data"]]
row.names(cell_type)<-row.names(adult_svz_caudate_umap_coordinates)

adult_neocortex_umap_coordinates <-adult_neocortex.umap[["data"]]
row.names(cell_type)<-row.names(adult_neocortex_umap_coordinates)

prenatal_gm_subclust_umap_coordinates <- prenatal_gm_subclust.umap[["data"]]
row.names(cell_type)<- row.names(prenatal_gm_subclust_umap_coordinates)

prenatal_cp_subclust_umap_coordinates <- prenatal_cp_subclust.umap[["layout"]]
prenatal_cp_subclust_umap_coordinates<-as.data.frame(prenatal_cp_subclust_umap_coordinates)
row.names(cell_type)<- row.names(prenatal_gm_subclust_umap_coordinates)

#merge UMAP coordinates and cell type annotations
combined <- cbind(prenatal_cp_all_umap_coordinates, cell_type)
combined <- cbind(adult_svz_caudate_umap_coordinates, cell_type)
combined <- cbind(adult_neocortex_umap_coordinates, cell_type)
combined <- cbind(prenatal_gm_subclust_umap_coordinates, cell_type)
combined <- cbind(prenatal_cp_subclust_umap_coordinates, cell_type)

colnames(combined)[1]<-"UMAP1"
colnames(combined)[2]<-"UMAP2"

#function for plotting coordinates of umap plot where input umap_df is the dataframe
#holding UMAP x coordinates in column UMAP1 and umap Y coordinates in UMAP2
#plottitle is plot title
#color_var is variable for color ie. cell_type_name (cell type label) or FABP7 expression
plot_umap <-function (umap_df, plottitle){
   p <- umap_df %>%
    ggplot(aes(x = UMAP1, 
               y = UMAP2, 
               color = cell_type_name
    ))+
    geom_point()+
    labs(x = "UMAP1",
         y = "UMAP2",
         subtitle = paste(plottitle))
   plot(p)
   
  #ggsave(paste(plottitle, ".png"))
   #save as png file, to save as pdf write "pdf", or jpg, write "jpg"; res = resolution
   png(filename=paste(plottitle,".png"), width=2650, height=2000, res=300)
   #function that makes the plot
   plot(p)
   dev.off()
}

plot_umap_FABP7 <-function (umap_df, plottitle){
  p <- umap_df %>%
    ggplot(aes(x = UMAP1, 
               y = UMAP2, 
               color = FABP7
    ))+
    geom_point()+
    labs(x = "UMAP1",
         y = "UMAP2",
         subtitle = paste(plottitle))
  plot(p)
  
  #ggsave(paste(plottitle, ".png"))
  #save as png file, to save as pdf write "pdf", or jpg, write "jpg"; res = resolution
  png(filename=paste(plottitle,".png"), width=2650, height=2000, res=300)
  #function that makes the plot
  plot(p)
  dev.off()
}

plot_umap(umap_df = combined, plottitle = "Ramos et al. (2022) All Prenatal Cortical Plate UMAP")
plot_umap(umap_df = combined, plottitle = "Ramos et al. (2022) Adult SVZ+caudate samples UMAP")
plot_umap(umap_df = combined, plottitle = "Ramos et al. (2022) Adult Neocortex samples UMAP")
plot_umap(umap_df = combined, plottitle = "Ramos et al. (2022) Prenatal Germinal Matrix Subclustered UMAP")
plot_umap(umap_df = combined, plottitle = "Ramos et al. (2022) Prenatal Cortical Plate Subclustered UMAP")

#function for extracting FABP7 expression from the RNA seq dataframes
#and merging with the cell types dataframe
#input rna_seq_df is the RNA seq dataframe ie. "prenatal_gm_all_full" dataframe
#input cell_type_df is the cell type dataframe containing the cell types
fabp7_processing<-function(rna_seq_df, cell_type_df){
  FABP7<-as.data.frame(rna_seq_df["FABP7",])
  FABP7_t <- t(FABP7)
  row.names(cell_type_df)<- colnames(rna_seq_df)
  combined <- cbind(FABP7_t, cell_type_df)
  return(combined)
}

#extract FABP7 expression and combine cell type annotations
combined_prenatal_cp_all <-fabp7_processing(rna_seq_df = prenatal_cp_all_full, cell_type_df = cell_type_prenatal_cp_all)
combined_prenatal_cp_subclustered <- fabp7_processing(rna_seq_df = prenatal_cp_subclust_full, cell_type_df = cell_type_prenatal_cp_subclust)
combined_prenatal_gm_all <- fabp7_processing(rna_seq_df = prenatal_gm_all_full, cell_type_df = cell_type_prenatal_gm_all)
combined_prenatal_gm_subclust<- fabp7_processing(rna_seq_df = prenatal_gm_subclust_full, cell_type_df = cell_type_prenatal_gm_subclustered)
combined_adult_neocortex<-fabp7_processing(rna_seq_df = adult_neocortex_full, cell_type_df = cell_type_adult_neocortex)
combined_adult_svz_caudate <- fabp7_processing(rna_seq_df = adult_svz_caudate_full, cell_type_df = cell_type_adult_svz_caudate)

#function for plotting FABP7 histogram from the different Ramos et al scRNA seq matrices
#input combined_fabp7 is the dataframe extracted for FABP7 expression and which also contains a column
#for cell type
plot_fabp7_histogram<-function(combined_fabp7, plottitle){
  p <-ggplot(combined_fabp7, aes(x = cell_type_name, 
                                           y = FABP7, color = cell_type_name))+
    geom_boxplot()+
    geom_jitter(color = "black", size = 0.7, alpha=0.9) +
    labs(x = "Cell Type",
         y = "FABP7 Expression (Combined Log-Normalized Average Expression)",
         subtitle = paste(plottitle))
  
  plot(p) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

#plot FABP7 
plot_fabp7_histogram(combined_fabp7 = combined_prenatal_gm_all, plottitle = "Ramos et al (2022) FABP7 Expression in Prenatal Germinal Matrix All")
plot_fabp7_histogram(combined_fabp7 = combined_prenatal_gm_subclust, plottitle = "Ramos et al (2022) FABP7 Expression in Prenatal Germinal Matrix Subclustered")
plot_fabp7_histogram(combined_fabp7 = combined_prenatal_cp_subclustered, plottitle = "Ramos et al (2022) FABP7 Expression in Prenatal Cortical Plate Subclustered")
plot_fabp7_histogram(combined_fabp7 = combined_prenatal_cp_all, plottitle = "Ramos et al (2022) FABP7 Expression in Prenatal Cortical Plate All")
plot_fabp7_histogram(combined_fabp7 = combined_adult_neocortex, plottitle = "Ramos et al (2022) FABP7 Expression in Adult Neocortex")
plot_fabp7_histogram(combined_fabp7 = combined_adult_svz_caudate, plottitle = "Ramos et al (2022) FABP7 Expression in Adult Subventricular Zone + Caudate")

 
p <-ggplot(combined_prenatal_gm_all, aes(x = cell_type_name, 
                           y = FABP7, color = cell_type_name))+
  geom_boxplot()+
  geom_jitter(color = "black", size = 0.7, alpha=0.9) +
  labs(x = "Cell Type",
       y = "FABP7 Expression (Combined Log-Normalized Average Expression",
       subtitle = paste("Ramos et al (2022) FABP7 Expression in Prenatal Germinal Matrix All"))

plot(p) + theme(axis.text.x = element_text(angle = 45, hjust = 1))


#plot FABP7 expression over the UMAP plot
FABP7_prenatal_cp <- as.data.frame(prenatal_cp_all_full["FABP7",])
FABP7_prenatal_cp_t <- t(FABP7_prenatal_cp)

FABP7_prenatal_cp_all <- cbind(combined, FABP7_prenatal_cp_t)

plot_umap(umap_df = FABP7_prenatal_cp_all, plottitle = "Ramos et al. (2022) All Prenatal Cortical Plate UMAP - FABP7 Expression")

FABP7_adult_svz_caudate <- adult_svz_caudate["FABP7",]
FABP7_adult_svz_caudate_t <- t(FABP7_adult_svz_caudate)
FABP7_adult_sv_caudate <- cbind(combined, FABP7_adult_svz_caudate_t)
plot_umap_FABP7(umap_df = FABP7_adult_sv_caudate, plottitle = "Ramos et al. (2022) Adult SVZ + Caudate UMAP - FABP7 Expression")

FABP7_adult_neocortex <-adult_neocortex["FABP7",]
FABP7_adult_neocortex_t <-t(FABP7_adult_neocortex)
FABP7_adult_neocortex <- cbind(combined, FABP7_adult_neocortex_t)
plot_umap_FABP7(umap_df = FABP7_adult_neocortex, plottitle = "Ramos et al. (2022) Adult Neocortex UMAP - FABP7 Expression")

FABP7_prenatal_gm_subclust <-prenatal_gm_subclust_full["FABP7",]
FABP7_prenatal_gm_subclust_t <-t(FABP7_prenatal_gm_subclust)
FABP7_prenatal_gm_subclust <- cbind(combined, FABP7_prenatal_gm_subclust_t)
plot_umap_FABP7(umap_df = FABP7_prenatal_gm_subclust, plottitle = "Ramos et al. (2022) Prenatal Germinal Matrix Subclustered - FABP7 Expression")

FABP7_prenatal_cp_subclust <- prenatal_cp_subclust["FABP7",]
FABP7<- t(FABP7_prenatal_cp_subclust)
FABP7_prenatal_cp_subclust<- cbind(combined, FABP7)
plot_umap_FABP7(umap_df = FABP7_prenatal_cp_subclust, plottitle = "Ramos et al. (2022) Prenatal Cortical Plate Subclustered - FABP7 Expression")

#run this command if error "Error in app$vspace(new_style$`margin-top` %||% 0)
#: attempt to apply non-function" pops up
install.packages("cli")
install.packages("palmerpenguins")

library(tidyverse)
library(palmerpenguins)

#https://datavizpyr.com/how-to-make-umap-plot-in-r/

penguins <- penguins %>% 
  drop_na() %>%
  select(-year)%>%
  mutate(ID=row_number()) 

penguins_meta <- penguins %>%
  select(ID, species, island, sex)


############## try with the penguin method

prenatal_gm_t<- t(prenatal_gm_all_full)

prenatal_gm_t<-as.data.frame(prenatal_gm_t)

prenatal_gm_t$sample <- row.names(prenatal_gm_t)
prenatal_gm_t$celltype<-cell_type$cell_type_name

prenatal_gm_t <- prenatal_gm_t %>%
  drop_na() %>%
  mutate(ID = row_number())

prenatal_gm_meta <- prenatal_gm_t %>%
  select(ID, sample, celltype)

row.names(prenatal_gm_meta)<-NULL

#perform UMAP

row.names(prenatal_gm_t) <- NULL

set.seed(142)

# umap_prenatal_gm <-prenatal_gm_t%>%
#   select(where(is.numeric)) %>%
#   column_to_rownames("ID") %>%
#   scale() %>%
#   umap()
#   





prenatal_gm_t <- prenatal_gm_t[,1:27612]
prenatal_gm_t_scaled <- scale(prenatal_gm_t)

prenatal_gm_all_combined <- list(expr = prenatal_gm_t, cellnames = cell_type$cell_type_name)

prenatal_gm_umap <-umap(prenatal_gm_all_combined$expr, labels=as.factor(prenatal_gm_all_combined$cellnames), controlscale=TRUE, scale=3)

umap_df <- prenatal_gm_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(prenatal_gm_meta, by="ID")

umap_df %>% head()

#         UMAP1      UMAP2 ID         sample                       celltype
# 1 -0.72393588  0.3250972  1       RNA.0.UD Undefined cell type or lineage
# 2  2.27413533 -0.2562321  2       RNA.1.IN                    Interneuron
# 3  0.06718826  0.1993333  3 RNA.2.L2.3.CPN     Cortical projection neuron
# 4  2.30525370 -0.8130698  4       RNA.3.IN                    Interneuron
# 5  2.68479496 -0.3473288  5       RNA.4.IN                    Interneuron
# 6  1.69470744 -0.5814841  6       RNA.5.IN                    Interneuron

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = celltype
             ))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "Prenatal GM All UMAP plot")
ggsave("UMAP_plot_prenatal_gm_all.png")

FABP7_expr <-prenatal_gm_t[,"FABP7"]

FABP7_expr<- as.data.frame(FABP7_expr)

combined_umap_df<-cbind(umap_df, FABP7_expr)

combined_umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = FABP7_expr
  ))+
    geom_point()+
    labs(x = "UMAP1",
         y = "UMAP2",
         subtitle = "FABP7 Expression in Ramos Prenatal GM All UMAP plot")

