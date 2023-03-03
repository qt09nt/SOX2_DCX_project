#Ramos et al (2022) glioblastoma dataset obtained from the supplementary information section of
#"An atlas of late prenatal human neurodevelopment resolved by single-nucleus transcriptomics"
#https://www.nature.com/articles/s41467-022-34975-2

#Seurat tutorial at 
#satijalab.org/seurat/articles/pbmc3k_tutorial.html

#In-depth-NGS-Data-Analysis-Course
#https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_marker_identification.html

setwd("C:\\Users\\qt09n\\Desktop\\Technical Analyst I UHN May 4 2021\\organoid group\\Sofia\\late prenatal human neurodevelopment resolved by single nucleus transcriptomics\\output")

#GSM672082 sample meta info:
#Sample 1,	Brain, cortical plate,	17 weeks gestation,	 female

#this file is results of tweaking Seurat tutorial settings with Alberto's suggested parameter changes on March 1st
GSM6720852<-readRDS("GSM6720852_march_1_2023_settings.rds")

#this file is all the pre-processing done with default Seurat tutorial settings
#GSM6720852<-readRDS(file = "GSM6720852.rds")

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

#March 1 
#try changing the filtering step parameters
GSM6720852 <- subset(GSM6720852, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)


########### Normalizing the data

#Normalized values are stored in GSM6720852[["RNA"]]@data
GSM6720852 <- NormalizeData(GSM6720852, normalization.method = "LogNormalize", scale.factor = 10000)

GSM6720852 <- NormalizeData(GSM6720852)

################### Identification of highly variable features (feature selection)

#we next calculate a subset of features that exhibit high cell-to-cell variation (ie. they are highly expressed
#in some cells, and lowly expressed in others).
#FindVariableFeaturs() function; by default returns 2000 features per dataset which will be used in downstream
#analysis ie. PCA

GSM6720852<- FindVariableFeatures(GSM6720852, selection.method = "vst", nfeatures = 2000)

#March 1 2023
#change the FindVariableFeatures parameters
GSM6720852<- FindVariableFeatures(GSM6720852, selection.method = "vst", nfeatures = 5000)


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

#FindClutsers() function implements process of applying modularity optimization techniques i.e 
#Louvain algorithm (default) or SLM (SLM, Blondel et al., Journal of statistical Mechanics),
#to iteratively group cells together with goal of optimizing standard modularity function


#the resolution parameter sets the 'granularity' of the downstream clustering, with increased values
#leading to a greater number of clusters
#usually a setting of this parameter between 0.4 to 1.2 typically returns good results for single-cell 
#datasets of around 3K cells. Optimal resolution often increases for larger datasets. 


GSM6720852 <- FindClusters(GSM6720852, resolution = 0.8)

#March 1 try out different values for resolution parameter
GSM6720852 <- FindClusters(GSM6720852, resolution = 1.0)

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

######## March 1 
# AAACCCAAGCATCAAA-1 AAACCCAAGTTTCGAC-1 AAACCCACACAAACGG-1 AAACCCACACCTGTCT-1 AAACCCACATACACCA-1 
#                 0                 15                 11                 14                  0 
# Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20

############# Run non-linear dimensional reduction (UMAP/tSNE)

GSM6720852 <- RunUMAP(GSM6720852, dims = 1:10)

#note that you can set 'label = TRUE' or use the LabelClusters function to help label
#individual clusters
DimPlot(GSM6720852, reduction = "umap", label.size = 5, label = TRUE)

setwd("C:\\Users\\qt09n\\Desktop\\Technical Analyst I UHN May 4 2021\\organoid group\\Sofia\\late prenatal human neurodevelopment resolved by single nucleus transcriptomics\\output")
saveRDS(GSM6720852, file = "GSM6720852.rds")
saveRDS(GSM6720852, file = "GSM6720852_march_1_2023_settings.rds")

#finding differentially expressed features (cluster biomarkers)

#find all markers of cluster 10
cluster10.markers <- FindMarkers(GSM6720852, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 5)

# |+++++++                                           | 14% ~09m 26s      Error in paste(deparse(expr, width.cutoff, ...), collapse = collapse) : 
#   could not allocate memory (1 Mb) in C function 'R_AllocStringBuffer'

#find all markers distinguishing cluster 10 from cluster 11 and cluster 14
cluster10.markers <- FindMarkers(GSM6720852, ident.1 = 10, ident.2 = c(11, 14), min.pct = 0.25)
head(cluster10.markers, n = 5)

# The results table output contains the following columns:
#   
# p_val: p-value not adjusted for multiple test correction
# avg_logFC: average log2 fold change. Positive values indicate that the gene is more highly expressed in the cluster.
# pct.1: The percentage of cells where the gene is detected in the cluster
# pct.2: The percentage of cells where the gene is detected on average in the other clusters
# p_val_adj: Adjusted p-value, based on bonferroni correction using all genes in the dataset, used to determine significance
# cluster: identity of cluster
# gene: Ensembl gene ID
# symbol: gene symbol
# biotype: type of gene
# description: gene description


#                    p_val avg_log2FC pct.1 pct.2     p_val_adj
# GRAMD1B    1.876981e-184   2.992480 0.955 0.465 4.704842e-180
# AL589740.1 2.433839e-172   2.490348 0.997 0.888 6.100660e-168
# MEIS2      1.426553e-163   1.842773 0.991 0.829 3.575799e-159
# NAV3       6.207860e-161  -2.766396 0.599 0.937 1.556062e-156
# TENM4      1.709841e-158   2.457994 0.949 0.579 4.285887e-154

#March 1 results 

#          p_val avg_log2FC pct.1 pct.2 p_val_adj
# IGSF21       0  1.9911392 0.832 0.052         0
# PTPRU        0  0.6658501 0.357 0.007         0
# GRIK3        0  2.4019151 0.968 0.102         0
# OLFM3        0  1.6455130 0.730 0.056         0
# TMEM178A     0  1.7522125 0.867 0.089         0


#find markers for every cluster compared to all remaining cells, report only the positive ones
GSM6720852.markers <- FindAllMarkers(GSM6720852, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(GSM6720852.markers, "GSM6720852_markers_march_1_settings.rds")

GSM6720852.markers <- FindAllMarkers(GSM6720852, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 2)

GSM6720852.markers %>%
  group_by(cluster) %>%
  slice_max(n =2, order_by = avg_log2FC)

# A tibble: 14 x 7
# Groups:   cluster [7]
#       p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene      
# <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>     
#   1 0               2.90 0.938 0.082 0         9       CUX2      
# 2 0               2.68 0.996 0.56  0         9       RBFOX1    
# 3 0               2.56 0.955 0.215 0         10      GRAMD1B   
# 4 0               2.48 0.645 0.043 0         10      RBFOX3    
# 5 0               4.46 0.852 0.046 0         11      ADARB2    
# 6 0               3.64 0.942 0.128 0         11      PDZRN3    
# 7 0               2.77 0.933 0.04  0         12      PRDM16    
# 8 9.23e-166       2.40 0.475 0.073 2.31e-161 12      TMEM132D  
# 9 0               4.03 0.953 0.079 0         13      NXPH1     
# 10 2.84e-191       3.48 1     0.484 7.12e-187 13      NRXN3     
# 11 0               3.86 0.69  0.066 0         14      HS3ST4    
# 12 0               2.54 0.624 0.035 0         14      CLSTN2    
# 13 6.04e- 53       3.76 0.468 0.051 1.51e- 48 15      AC007402.1
# 14 1.24e- 99       3.66 0.516 0.034 3.12e- 95 15      PCDH15  


### March 1 2023 settings

#        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene       
#       <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>      
# 1 0               4.11 1     0.608 0         1       NRXN3      
# 2 0               3.76 0.845 0.224 0         1       PDZRN3     
# 3 0               2.90 0.982 0.285 0         4       EEPD1      
# 4 0               2.66 0.931 0.167 0         4       TNC        
# 5 4.66e-186       2.15 0.831 0.254 1.17e-181 5       GPR39      
# 6 1.90e-178       2.04 0.993 0.661 4.76e-174 5       GRAMD1B    
# 7 8.33e-221       2.29 0.874 0.221 2.09e-216 8       AC016205.1 
# 8 6.35e-234       2.17 0.932 0.241 1.59e-229 8       DIAPH3     
# 9 0               4.52 1     0.122 0         10      HS3ST4     
# 10 1.16e-266       3.93 0.997 0.343 2.90e-262 10      DPP10      
# 11 1.58e-186       2.85 0.98  0.342 3.96e-182 11      HS6ST3     
# 12 6.51e-116       2.71 0.698 0.202 1.63e-111 11      IQCJ-SCHIP1
# 13 3.01e- 12       2.08 0.527 0.309 7.53e-  8 12      KCNIP4     
# 14 7.59e-224       3.18 0.957 0.147 1.90e-219 13      APOLD1     
# 15 1.92e-162       2.74 0.995 0.26  4.81e-158 13      DIAPH3     
# 16 4.81e-120       3.02 0.972 0.328 1.21e-115 15      SAMD4A     
# 17 9.91e-106       2.66 0.983 0.455 2.48e-101 15      SHROOM3    
# 18 2.04e-138       3.63 0.949 0.208 5.10e-134 16      KCNQ5      
# 19 1.54e- 85       3.14 0.942 0.401 3.87e- 81 16      RALYL      
# 20 8.90e- 44       3.07 1     0.565 2.23e- 39 17      KIF26B     
# 21 2.34e- 31       2.89 0.671 0.195 5.88e- 27 17      TMEM132D   
# 22 7.62e-148       4.61 0.849 0.072 1.91e-143 18      PCDH15     
# 23 1.38e- 58       4.45 0.658 0.108 3.46e- 54 18      AC007402.1 
# 24 4.23e- 55       4.62 0.784 0.09  1.06e- 50 19      GALNTL6    
# 25 1.00e- 56       3.26 0.811 0.088 2.52e- 52 19      GABRG3     
# 26 6.45e- 80       5.08 1     0.056 1.62e- 75 20      LGMN       
# 27 7.10e-214       5.08 1     0.018 1.78e-209 20      SPP1 
#Seurat can test for differential expression with the test.use parameter. ie. the ROC test returns the 
#'classification power' for any individual marker (ranging from 0 - random, and 1 - perfect)

#March 1st settings using
#GSM6720852.markers <- FindAllMarkers(GSM6720852, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Rearrange order of columns to make clearer
GSM6720852.markers<- GSM6720852.markers[,c(6,7,1:5)]

#Save all the identified markers to a file
write.csv(GSM6720852.markers, "GSM6720852.csv", quote = F)

#return top 10 markers for clusters 1:20 and create a dataframe
GSM6720852_top10_markers <-GSM6720852.markers %>%
  group_by(cluster) %>%
  slice_max(n =10, order_by = avg_log2FC) %>%    #10 gene markers per cluster
  print(n = 210)                                 #10 genes x 21 clusters = 210

#write these results to a file
write.csv(GSM6720852_top10_markers, "GSM6720852_top10_markers.csv", quote = F)

# A tibble: 42 x 7
# Groups:   cluster [21]
#         p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene   Cell Type
#         <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>      
# 1   0               1.90 0.994 0.692 0         0       SORBS2     
# 2   6.11e-272       1.89 0.721 0.197 1.53e-267 0       EML6       
# 3   0               4.11 1     0.608 0         1       NRXN3      
# 4   0               3.76 0.845 0.224 0         1       PDZRN3     
# 5   1.27e- 10       1.25 0.522 0.651 3.19e-  6 2       MARCKS     
# 6   3.90e-123       1.23 0.935 0.973 9.77e-119 2       MAP1B      
# 7   1.79e- 13       1.19 0.603 0.777 4.48e-  9 3       TUBA1A     
# 8   1.54e-  4       1.16 0.483 0.652 1   e+  0 3       MARCKS     
# 9   0               2.90 0.982 0.285 0         4       EEPD1      
# 10  0               2.66 0.931 0.167 0         4       TNC        
# 11  4.66e-186       2.15 0.831 0.254 1.17e-181 5       GPR39      
# 12  1.90e-178       2.04 0.993 0.661 4.76e-174 5       GRAMD1B    
# 13  4.65e-123       1.57 0.965 0.504 1.16e-118 6       DSCAM      
# 14  1.56e-116       1.57 0.807 0.334 3.92e-112 6       NRP1       
# 15  2.44e-150       1.97 0.902 0.338 6.12e-146 7       HS6ST3     
# 16  8.36e-151       1.85 0.995 0.565 2.09e-146 7       SATB2      
# 17  8.33e-221       2.29 0.874 0.221 2.09e-216 8       AC016205.1 
# 18  6.35e-234       2.17 0.932 0.241 1.59e-229 8       DIAPH3     
# 19  1.18e-127       1.81 0.969 0.851 2.96e-123 9       MT-CO1     
# 20  1.13e-120       1.79 0.944 0.823 2.84e-116 9       MT-CO2     
# 21  0               4.52 1     0.122 0         10      HS3ST4     
# 22  1.16e-266       3.93 0.997 0.343 2.90e-262 10      DPP10      
# 23  1.58e-186       2.85 0.98  0.342 3.96e-182 11      HS6ST3     
# 24  6.51e-116       2.71 0.698 0.202 1.63e-111 11      IQCJ-SCHIP1
# 25  3.01e- 12       2.08 0.527 0.309 7.53e-  8 12      KCNIP4     
# 26  2.59e- 21       1.67 0.393 0.153 6.49e- 17 12      CLSTN2     
# 27  7.59e-224       3.18 0.957 0.147 1.90e-219 13      APOLD1     
# 28  1.92e-162       2.74 0.995 0.26  4.81e-158 13      DIAPH3     
# 29  7.79e- 72       1.61 1     0.639 1.95e- 67 14      NRXN3      
# 30  1.75e- 64       1.50 1     0.762 4.37e- 60 14      ERBB4      
# 31  4.81e-120       3.02 0.972 0.328 1.21e-115 15      SAMD4A     
# 32  9.91e-106       2.66 0.983 0.455 2.48e-101 15      SHROOM3    
# 33  2.04e-138       3.63 0.949 0.208 5.10e-134 16      KCNQ5      
# 34  1.54e- 85       3.14 0.942 0.401 3.87e- 81 16      RALYL      
# 35  8.90e- 44       3.07 1     0.565 2.23e- 39 17      KIF26B     
# 36  2.34e- 31       2.89 0.671 0.195 5.88e- 27 17      TMEM132D   
# 37  7.62e-148       4.61 0.849 0.072 1.91e-143 18      PCDH15     
# 38  1.38e- 58       4.45 0.658 0.108 3.46e- 54 18      AC007402.1 
# 39  4.23e- 55       4.62 0.784 0.09  1.06e- 50 19      GALNTL6    
# 40  1.00e- 56       3.26 0.811 0.088 2.52e- 52 19      GABRG3     
# 41  6.45e- 80       5.08 1     0.056 1.62e- 75 20      LGMN       
# 42  7.10e-214       5.08 1     0.018 1.78e-209 20      SPP1  

View(GSM6720852.markers)

cluster0.markers <- FindMarkers(GSM6720852, ident.1 = 0, logfc.threshold = 2, test.use = "roc", only.pos = TRUE)

# Warning message:
#   In FindMarkers.default(object = data.use, slot = data.slot, counts = counts,  :
#                            No features pass logfc.threshold threshold; returning empty data.frame


#Visualizing marker Expression. VlnPlot() shows expression probability distributions 
#across clusters), and FeaturePlot() visualizes feature expression on a tSNE or PCA plot 
#are the most commonly used visualization.
#other suggestions: RidgePlot(), CellScatter(), DotPlot() as other ways to view dataset
VlnPlot(GSM6720852, features = c("FABP7", "SOX2"))
VlnPlot(GSM6720852, features = c("FABP7"))

#can plot raw counts as well

VlnPlot(GSM6720852, features = c("FABP7", "DCX"), slot = "counts", log = TRUE)
VlnPlot(GSM6720852, features = c("FABP7"), slot = "counts", log = TRUE)

FeaturePlot(GSM6720852, features = c("FABP7"))

#DoHeatmap() generates an expression heatmap for given cells and features
#plot the top 20 markers (or all markers if less than 20) for each cluster

#hard to read; too many genes
GSM6720852.markers %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)-> top10
DoHeatmap(GSM6720852, features = top10$gene) +NoLegend()

GSM6720852.markers %>% 
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)-> top5
DoHeatmap(GSM6720852, features = top5$gene) +NoLegend()


#doesn't work:
DimPlot(
  GSM6720852,
  "tnse",
  label.size = 8) +
  ggtitle("tSNE")

####################### https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_marker_identification.html

######### Get an idea of cell type identity by expression of different identified markers by
#cluster using the FeaturePlot() function

FeaturePlot(GSM6720852, features=c("FABP7"))

#Assigning Cell Type Identity to Clusters


