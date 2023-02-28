#Ramos et al (2022) glioblastoma dataset obtained from the supplementary information section of
#"An atlas of late prenatal human neurodevelopment resolved by single-nucleus transcriptomics"
#https://www.nature.com/articles/s41467-022-34975-2

#Seurat tutorial at 
#satijalab.org/seurat/articles/pbmc3k_tutorial.html


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

#Normalized values are stored in GSM6720852[["RNA"]]@data
GSM6720852 <- NormalizeData(GSM6720852, normalization.method = "LogNormalize", scale.factor = 10000)

GSM6720852 <- NormalizeData(GSM6720852)

################### Identification of highly variable features (feature selection)

#we next calculate a subset of features that exhibit high cell-to-cell variation (ie. they are highly expressed
#in some cells, and lowly expressed in others).

GSM6720852<- FindVariableFeatures(GSM6720852, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(GSM6720852), 10)

top10
# [1] "ADARB2"     "HS3ST4"     "AC007402.1" "PDZRN3"     "NRXN3"      "NXPH1"      "HBA1"       "TMEM132D"  
# [9] "HBA2"       "KCNIP4"

#plot variable features with and without labels
plot1 <- VariableFeaturePlot(GSM6720852)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2 

#need to expand the Plots viewport in R Studio to prevent the error below:
# Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : 
#                      Viewport has zero dimension(s)
#                    In addition: Warning messages:
#                      1: Transformation introduced infinite values in continuous x-axis 
#                    2: Removed 127 rows containing missing values (geom_point). 
#                    3: Transformation introduced infinite values in continuous x-axis 
#                    4: Removed 127 rows containing missing values (geom_point). 

####### Scaling the data 
all.genes <- rownames(GSM6720852)
GSM6720852 <- ScaleData(GSM6720852, features = all.genes)
GSM6720852 <- ScaleData(GSM6720852)

GSM6720852 <- ScaleData(GSM6720852, vars.to.regress = "percent.mt")



############### Perform linear dimensional reduction

GSM6720852 <- RunPCA(GSM6720852, features = VariableFeatures(object = GSM6720852))

# PC_ 1 
# Positive:  RBFOX3, CUX2, DSCAML1, DIP2C, CDH4, PLXNA2, PHF21B, AC092683.1, GTF2IRD1, SH3PXD2A 
# KDM4B, CORO2B, MIAT, LINC00854, PTPRN2, PTPRS, HIVEP3, CTIF, LDLRAD4, RAB11FIP3 
# CASTOR2, SLC22A23, TNRC18, NFASC, SHANK2, ARL17B, SH3RF3, NXN, CACNA1B, AGAP1 
# Negative:  ERBB4, AL589740.1, SPP1, TNC, SLC7A11, HBA2, LRP1B, BCAS1, LYZ, HLA-DRB1 
# AC074351.1, COL1A1, PPP1R3G, CD8B2, AC087636.1, MGAT4C, LINC00844, AREG, RNASE1, DMBT1 
# AC084757.4, DMKN, CAVIN3, APOBEC1, AL121760.1, SLC18A1, AL390061.1, MS4A14, SCGB1D2, AL358075.2 
# PC_ 2 
# Positive:  ZNF423, PRDM16, CDH4, RBFOX3, KCNN3, TTYH2, GRAMD1B, PREX1, GPR39, AL589740.1 
# TENM4, SLC22A23, NEK6, CLMP, CA12, CCDC88C, KCNQ3, KIF26B, FGFR2, NRG1 
# GLI2, LRP8, MIR9-3HG, NFIC, EPHB1, SLC39A11, PTK7, GNG4, SH3PXD2A, MIRLET7BHG 
# Negative:  ST8SIA5, ADARB2, ZNF536, ZBTB16, NPAS1, RIPOR2, CACNA1C, CELF4, SDK2, NRXN3 
# DLX6-AS1, PDZRN3, GAD2, FGD3, KIRREL3, THRB, DSCAML1, GAD1, ERBB4, KCNIP1 
# GRIP2, SLC6A1, PRKCA, PLS3, BRINP2, SORCS3, NR2F2-AS1, SLIT1, TNR, ATRNL1 
# PC_ 3 
# Positive:  CUX2, DSCAM, RBFOX1, KLHL29, NFASC, CACNA1E, MPPED1, RBFOX3, SATB2, NRP1 
# SLA, PLXNA4, FRMD4B, CLMP, DNAH10, DCC, HS6ST3, MYT1L, ATP8A2, PLXNA2 
# MCTP1, CACNA1A, NELL2, CACNB4, KCNK2, PTPRD, SEMA6D, SPHKAP, DOPEY2, SNTG2 
# Negative:  PRDM16, CELSR1, MMD2, COL27A1, GLI2, NEK6, BOC, TBL1X, FAM182B, SLC6A1 
# TTYH1, BMP7, LIPG, AC087564.1, ABTB2, DAPK2, ST5, MAGI1, LRIG1, SEMA5B 
# DOCK1, EGFR, NPAS3, LTBP1, LDLR, CDH23, SYT17, PDGFRB, EEPD1, MIR9-3HG 
# PC_ 4 
# Positive:  ZNF536, ST8SIA5, ADARB2, ELMO1, TTYH2, FGD3, NPAS1, GNG4, THRB, PDZRN3 
# DLX6-AS1, ZNF423, ZBTB16, GAD2, SMOC1, RIPOR2, CORO2B, SERINC5, RBFOX3, PHF21B 
# LRP8, KCNH8, CCDC88C, DAPK1, GRAMD1B, SH3PXD2A, EML6, SIPA1L2, GSE1, AC106869.1 
# Negative:  NTRK3, KAZN, HS3ST4, DAB1, MYO5B, CSMD1, CALN1, PRKCE, KCNMA1, SLIT3 
# ASTN2, FRMPD4, CCBE1, DOCK9, PRICKLE2, RGS6, HS6ST3, GAS7, MTUS2, LIMCH1 
# CLSTN2, VLDLR-AS1, SYN2, FAT3, MCTP1, FAM160A1, CSGALNACT1, SPHKAP, ITPR1, KIAA1217 
# PC_ 5 
# Positive:  FAM135B, NXPH1, PTPRT, TENM3, BRINP2, NXPH2, GRIK3, LHFPL3, ELFN1, LHX6 
# MAF, TOX2, TENM2, EPHA6, GALNT17, GRIA4, KCNC2, PTCHD4, NELL1, GAD1 
# SLC4A4, ARX, PRKD1, TENM1, RUNX1T1, SYT16, CNTNAP2, CNTNAP4, CUX1, FAM189A1 
# Negative:  ADARB2, SORCS3, NR2F2-AS1, KIRREL3, TNR, SLIT1, ZBTB16, THRB, NR2F2, NR3C2 
# CALB2, KCNJ3, SDK2, PRKCA, RXRA, KITLG, PDZRN3, DPF3, PIP5K1B, EYA1 
# PROX1, ST8SIA5, SHISA9, SH3RF3, GLRA2, SPOCK1, U91319.1, PCDH9, SHANK2, LINC02398

print(GSM6720852[["pca"]], dims =1:5, nfeatures = 5)

# PC_ 1 
# Positive:  RBFOX3, CUX2, DSCAML1, DIP2C, CDH4 
# Negative:  ERBB4, AL589740.1, SPP1, TNC, SLC7A11 
# PC_ 2 
# Positive:  ZNF423, PRDM16, CDH4, RBFOX3, KCNN3 
# Negative:  ST8SIA5, ADARB2, ZNF536, ZBTB16, NPAS1 
# PC_ 3 
# Positive:  CUX2, DSCAM, RBFOX1, KLHL29, NFASC 
# Negative:  PRDM16, CELSR1, MMD2, COL27A1, GLI2 
# PC_ 4 
# Positive:  ZNF536, ST8SIA5, ADARB2, ELMO1, TTYH2 
# Negative:  NTRK3, KAZN, HS3ST4, DAB1, MYO5B 
# PC_ 5 
# Positive:  FAM135B, NXPH1, PTPRT, TENM3, BRINP2 
# Negative:  ADARB2, SORCS3, NR2F2-AS1, KIRREL3, TNR 

VizDimLoadings(GSM6720852, dims = 1:2, reduction = "pca")

DimPlot(GSM6720852, reduction="pca")

DimHeatmap(GSM6720852, dims =1, cells = 500, balanced = TRUE)
DimHeatmap(GSM6720852, dims =1:15, cells = 500, balanced = TRUE)

####### Determine the Dimensionality of the data

#Note: computationally expensive step, make take a long time; 
#can comment out for expediency
GSM6720852<-JackStraw(GSM6720852, num.replicate = 1000)
GSM6720852 <- ScoreJackStraw(GSM6720852, dims = 1:20)

ElbowPlot(GSM6720852)

############## Cluster the Cells

GSM6720852 <- FindNeighbors(GSM6720852, dims = 1:10)
GSM6720852 <- FindClusters(GSM6720852, resolution = 0.8)

# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
# 
# Number of nodes: 50854
# Number of edges: 1017130
# 
# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#    Maximum modularity in 10 random starts: 0.6922
#    Number of communities: 16
#    Elapsed time: 26 seconds

#Look at cluster IDs of the first 5 cells
head(Idents(GSM6720852), 5)

# AAACCCAAGAATTGCA-1 AAACCCAAGACAGCGT-1 AAACCCAAGACCCTTA-1 AAACCCAAGAGGCGGA-1 AAACCCAAGATGACCG-1 
#                 4                  4                  3                  2                  1 
# Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15

############# Run non-linear dimensional reduction (UMAP/tSNE)


