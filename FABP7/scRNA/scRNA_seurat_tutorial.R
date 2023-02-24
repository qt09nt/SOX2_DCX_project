#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Seurat - Guided Clustering Tutorial
# Compiled: January 11, 2022
# Source: vignettes/pbmc3k_tutorial.Rmd

# Setup the Seurat Object
# 
# For this tutorial, we will be analyzing the a dataset of Peripheral Blood Mononuclear Cells (PBMC) 
# freely available from 10X Genomics. There are 2,700 single cells that were sequenced on the Illumina NextSeq 500. 
# The raw data can be found here.
# # 
# # We start by reading in the data. The Read10X() function reads in the output of the cellranger pipeline from 10X, 
# returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules 
# for each feature (i.e. gene; row) that are detected in each cell (column).
# # 
# # We next use the count matrix to create a Seurat object. The object serves as a container that contains both data 
# (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. For a technical 
# discussion of the Seurat object structure, check out our GitHub Wiki. For example, the count matrix is stored 
# in pbmc[["RNA"]]@counts.

setwd("/Users/queenietsang/Desktop/ramos_scRNA/data/filtered_gene_bc_matrices/hg19")

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/queenietsang/Desktop/ramos_scRNA/data/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#What does data in a count matrix look like?
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
# 3 x 30 sparse Matrix of class "dgCMatrix"
# [[ suppressing 30 column names ‘AAACATACAACCAC-1’, ‘AAACATTGAGCTAC-1’, ‘AAACATTGATCAGC-1’ ... ]]

# CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
# TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
# MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .

# The . values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq 
# matrix are 0, Seurat uses a sparse-matrix representation whenever possible. This results in significant
# memory and speed savings for Drop-seq/inDrop/10x data.

dense.size <- object.size(as.matrix(pbmc.data))
dense.size
#709591472 bytes


sparse.size <- object.size(pbmc.data)
sparse.size
#29905192 bytes

dense.size/sparse.size
#23.7 bytes

# Standard pre-processing workflow
# 
# The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. 
#These represent the selection and filtration of cells based on QC metrics, data normalization 
#and scaling, and the detection of highly variable features.


#QC and selecting cells for further analysis

# Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria.
#A few QC metrics commonly used by the community include
# 
# The number of unique genes detected in each cell.
# Low-quality cells or empty droplets will often have very few genes
# Cell doublets or multiplets may exhibit an aberrantly high gene count
# Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
# The percentage of reads that map to the mitochondrial genome
# Low-quality / dying cells often exhibit extensive mitochondrial contamination
# We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts 
#originating from a set of features
# We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Where are QC metrics stored in Seurat?
#   
# The number of unique genes and total molecules are automatically calculated during CreateSeuratObject()
# You can find them stored in the object meta data

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

#                   orig.ident nCount_RNA nFeature_RNA percent.mt
# AAACATACAACCAC-1     pbmc3k       2419          779  3.0177759
# AAACATTGAGCTAC-1     pbmc3k       4903         1352  3.7935958
# AAACATTGATCAGC-1     pbmc3k       3147         1129  0.8897363
# AAACCGTGCTTCCG-1     pbmc3k       2639          960  1.7430845
# AAACCGTGTATGCG-1     pbmc3k        980          521  1.2244898

# In the example below, we visualize QC metrics, and use these to filter cells.
# 
# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalizing the data
# 
# After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a
# global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell 
# by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# Normalized values are stored in pbmc[["RNA"]]@data.

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# For clarity, in this previous line of code (and in future commands), we provide the default values for certain parameters
# in the function call. However, this isn’t required and the same behavior can be achieved with:

pbmc <- NormalizeData(pbmc)

# Identification of highly variable features (feature selection)
# 
# We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly 
# expressed in some cells, and lowly expressed in others). We and others have found that focusing on these genes in 
# downstream analysis helps to highlight biological signal in single-cell datasets.
# 
# Our procedure in Seurat is described in detail here, and improves on previous versions by directly modeling the 
# mean-variance relationship inherent in single-cell data, and is implemented in the FindVariableFeatures() function.
# By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data
# 
# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior 
#to dimensional reduction techniques like PCA. The ScaleData() function:
#   
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# This step takes too long! Can I make it faster?
#   
# Scaling is an essential step in the Seurat workflow, but only on genes that will be used as input to PCA. Therefore, 
# the default in ScaleData() is only to perform scaling on the previously identified variable features (2,000 by default). 
# To do this, omit the features argument in the previous function call, i.e.

pbmc <- ScaleData(pbmc)

# Your PCA and clustering results will be unaffected. However, Seurat heatmaps (produced as shown below with DoHeatmap()) 
# require genes in the heatmap to be scaled, to make sure highly-expressed genes don’t dominate the heatmap. To make sure we
# don’t leave any genes out of the heatmap later, we are scaling all genes in this tutorial. 
# 
# How can I remove unwanted sources of variation, as in Seurat v2?
#   
# In Seurat v2 we also use the ScaleData() function to remove unwanted sources of variation from a single-cell dataset. 
# For example, we could ‘regress out’ heterogeneity associated with (for example) cell cycle stage, or mitochondrial contamination.
# These features are still supported in ScaleData() in Seurat v3, i.e.:

pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

# However, particularly for advanced users who would like to use this functionality, we strongly recommend the use of our new 
# normalization workflow, SCTransform(). The method is described in our paper, with a separate vignette using Seurat v3 here.
# As with ScaleData(), the function SCTransform() also includes a vars.to.regress parameter. 

# Perform linear dimensional reduction
# 
# Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, 
# but can be defined using features argument if you wish to choose a different subset.

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# PC_ 1 
# Positive:  CST3, TYROBP, LST1, AIF1, FTL 
# Negative:  MALAT1, LTB, IL32, IL7R, CD2 
# PC_ 2 
# Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1 
# Negative:  NKG7, PRF1, CST7, GZMA, GZMB 
# PC_ 3 
# Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPA1 
# Negative:  PPBP, PF4, SDPR, SPARC, GNG11 
# PC_ 4 
# Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1 
# Negative:  VIM, IL7R, S100A6, S100A8, IL32 
# PC_ 5 
# Positive:  GZMB, FGFBP2, S100A8, NKG7, GNLY 
# Negative:  LTB, IL7R, CKB, MS4A7, RP11-290F20.3 

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")



# In particular DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and 
# can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and features
# are ordered according to their PCA scores. Setting cells to a number plots the ‘extreme’ cells on both ends of 
# the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find
# this to be a valuable tool for exploring correlated feature sets.

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
# 
# To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells
# based on their PCA scores, with each PC essentially representing a ‘metafeature’ that combines information across a 
# correlated feature set. The top principal components therefore represent a robust compression of the dataset. However, 
# how many components should we choose to include? 10? 20? 100?
# #   
# #   In Macosko et al, we implemented a resampling test inspired by the JackStraw procedure. We randomly permute a 
#   subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat 
#     this procedure. We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

# The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for 
# each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features
# with low p-values (solid curve above the dashed line). In this case it appears that there is a sharp drop-off in 
# significance after the first 10-12 PCs.

JackStrawPlot(pbmc, dims = 1:15)

# An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on 
# the percentage of variance explained by each one (ElbowPlot() function). In this example, we can observe an ‘elbow’
# around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.

ElbowPlot(pbmc)

# Identifying the true dimensionality of a dataset – can be challenging/uncertain for the user. We therefore suggest these three approaches to consider. The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example. The second implements a statistical test based on a random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is commonly used, and can be calculated instantly. In this example, all three approaches yielded similar results, but we might have been justified in choosing anything between PC 7-12 as a cutoff.
# 
# We chose 10 here, but encourage users to consider the following:
#   
#   Dendritic cell and NK aficionados may recognize that genes strongly associated with PCs 12 and 13 define rare immune subsets (i.e. MZB1 is a marker for plasmacytoid DCs). However, these groups are so rare, they are difficult to distinguish from background noise for a dataset of this size without prior knowledge.
# We encourage users to repeat downstream analyses with a different number of PCs (10, 15, or even 50!). As you will observe, the results often do not differ dramatically.
# We advise users to err on the higher side when choosing this parameter. For example, performing downstream analyses with only 5 PCs does significantly and adversely affect results.



# Finding differentially expressed features (cluster biomarkers)
# 
# Seurat can help you find markers that define clusters via differential expression. By default, 
# it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells.
# FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, 
# or against all cells.
# # 
# # The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, 
# and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between 
# the two groups. You can set both of these to 0, but with a dramatic increase in time - since this will test a large
# number of features that are unlikely to be highly discriminatory. As another option to speed up these computations, 
# max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is 
# set to. While there is generally going to be a loss in power, the speed increases can be significant and the most 
# highly differentially expressed features will likely still rise to the top.

## find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

