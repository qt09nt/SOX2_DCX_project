#Ramos et al (2022) glioblastoma dataset obtained from the supplementary information section of
#"An atlas of late prenatal human neurodevelopment resolved by single-nucleus transcriptomics"
#https://www.nature.com/articles/s41467-022-34975-2

library(dplyr)
library(Seurat)
library(patchwork)

#load the scRNAseq glioblastoma dataset 
GSM6720852_9C<-Read10X(data.dir = "D:/ramos/GSM6720852")
GSM6720853_9G<-Read10X(data.dir = "D:/ramos/GSM6720853")

#initialize seurat object with raw (non-normalized data)
GSM6720852<-CreateSeuratObject(counts = GSM6720852_9C, project = "ramos", min.cells = 3, min.features = 200)
GSM6720852

# An object of class Seurat 
# 25066 features across 58897 samples within 1 assay 
# Active assay: RNA (25066 features, 0 variable features)

#what does data in a count matrix look like?
#let's examine a few genes in the first 100 cells
GSM6720852_9C[c("FABP7", "SOX2", "DCX"), 1:100]
# 3 x 100 sparse Matrix of class "dgCMatrix"
# [[ suppressing 100 column names 'AAACCCAAGAAACACT-1', 'AAACCCAAGAAACCAT-1', 'AAACCCAAGAAACCCA-1' ... ]]
# 
# FABP7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# SOX2  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# DCX   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# 
# FABP7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# SOX2  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# DCX   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

#The . values in th ematrix represent 0s (no molecules detected). Since most values in an scRNAseq matrix are 0, 
#Seurat uses a sparse-matrix representation whenever possible. this results in significant memory and speed savings 
#for Drop-seq/inDrop/10x data

dense.size <- object.size(as.matrix(GSM6720852_9C))
# Error in asMethod(object) : 
#   Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 102

dense.size

sparse.size<-object.size(GSM6720852_9C)
sparse.size
#1094855952 bytes

dense.size/sparse.size

###################Standard pre-processing workflow
#the steps encompass the standard pre-processing workflow for scRNA-seq data in Seurat. These represent
#the selection and filtration of cells based on QC metrics, data normalization and scaling, and detection of
#highly variable features


###############QC and selecting cells for further analysis\

#Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats

GSM6720852[["percent.mt"]]<- PercentageFeatureSet(GSM6720852, pattern = "^MT-")

########Where are QC metrics stored in Seurat?
#the number of unique genes and totla molecules are automatically calculated during the CreateSeuratObject()
#you can find them stored in the object meta data

#Show QC metrics for the first 5 cells
head(GSM6720852@meta.data, 5)

#                     orig.ident nCount_RNA nFeature_RNA percent.mt
# AAACCCAAGAATTGCA-1      ramos        472          409  2.5423729
# AAACCCAAGACAGCGT-1      ramos        546          447  3.4798535
# AAACCCAAGACCCTTA-1      ramos        509          428  0.0000000
# AAACCCAAGAGGCGGA-1      ramos        557          454  0.5385996
# AAACCCAAGATGACCG-1      ramos        554          449  1.6245487

#in the example below, we visualize QC metrics and use these to filter cells
#we filter cells that have a unique feature counts over 2500 or less than 200
#we filter cells that have >5% mitochondrial counts
#visualize QC metrics as violiln plot
VlnPlot(GSM6720852, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#FeatureScatter is typically used to visualize the feature-feature relationships but can
#be used for anything calculated by the object, ie. columns in object metadata, PC scores
#etc.

plot1<-FeatureScatter(GSM6720852, feature1 = "nCount_RNA", feature2="percent.mt")
plot2<-FeatureScatter(GSM6720852, feature1= "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

GSM6720852 <- subset(GSM6720852, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

########### Normalizing the data

GSM6720852 <- NormalizeData(GSM6720852, normalization.method = "LogNormalize", scale.factor = 10000)

GSM6720852 <- NormalizeData(GSM6720852)

################### Identification of highly variable features (feature selection)

#we next calculate a subset of features that exhibit high cell-to-cell variation (ie. they are highly expressed
#in some cells, and lowly expressed in others).

GSM6720852<- FindVariableFeatures(GSM6720852, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GSM6720852), 10)

#plot variable features with and without labels
plot1 <- VariableFeaturePlot(GSM6720852)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2 + theme(legend.position="none")

# Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : 
#                      Viewport has zero dimension(s)
#                    In addition: Warning messages:
#                      1: Transformation introduced infinite values in continuous x-axis 
#                    2: Removed 127 rows containing missing values (geom_point). 
#                    3: Transformation introduced infinite values in continuous x-axis 
#                    4: Removed 127 rows containing missing values (geom_point). 


