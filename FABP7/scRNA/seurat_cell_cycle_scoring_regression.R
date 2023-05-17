#Seurat Cell-Cycle Scoring and Regression

library(Seurat)

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "C:/Users/Diamandis Lab II/Documents/Queenie/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
                      as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create our Seurat object and complete the initalization steps
marrow <- CreateSeuratObject(counts = exp.mat)
marrow <- NormalizeData(marrow)
marrow <- FindVariableFeatures(marrow, selection.method = "vst")
marrow <- ScaleData(marrow, features = rownames(marrow))

# If we run a PCA on our object, using the variable genes we found in FindVariableFeatures() above, 
# we see that while most of the variance can be explained by lineage, PC8 and PC10 are split on cell-cycle
# genes including TOP2A and MKI67. We will attempt to regress this signal from the data, so that
# cell-cycle heterogeneity does not contribute to PCA or downstream analysis.

marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 6:10, nfeatures.print = 10)

# PC_ 6 
# Positive:  SELL, ARL6IP1, CCL9, CD34, ADGRL4, BPIFC, NUSAP1, FAM64A, CD244, C030034L19RIK 
# Negative:  LY6C2, AA467197, CYBB, MGST2, ITGB2, PF4, CD74, ATP1B1, GP1BB, TREM3 
# PC_ 7 
# Positive:  HDC, CPA3, PGLYRP1, MS4A3, NKG7, UBE2C, CCNB1, NUSAP1, PLK1, FUT8 
# Negative:  F13A1, LY86, CFP, IRF8, CSF1R, TIFAB, IFI209, CCR2, TNS4, MS4A6C 
# PC_ 8 
# Positive:  WFDC17, SLC35D3, ADGRL4, VLDLR, CD33, H2AFY, P2RY14, IFI206, CCL9, CD34 
# Negative:  NUSAP1, UBE2C, KIF23, PLK1, CENPF, FAM64A, CCNB1, H2AFX, ID2, CDC20 
# PC_ 9 
# Positive:  SLC2A6, HBA-A1, HBA-A2, IGHV8-7, FCER1G, F13A1, HBB-BS, PLD4, HBB-BT, IGFBP4 
# Negative:  IGKC, JCHAIN, LY6D, MZB1, CD74, IGLC2, FCRLA, IGKV4-50, IGHM, IGHV9-1 
# PC_ 10 
# Positive:  CTSW, XKRX, PRR5L, RORA, MBOAT4, A630014C17RIK, ZFP105, COL9A3, CLEC2I, TRAT1 
# Negative:  H2AFX, FAM64A, ZFP383, NUSAP1, CDC25B, CENPF, GBP10, TOP2A, GBP6, GFRA1 

DimHeatmap(marrow, dims = c(8, 10))

# Assign Cell-Cycle Scores
# 
# First, we assign each cell a score, based on its expression of G2/M and S phase markers. 
# These marker sets should be anticorrelated in their expression levels, and cells expressing neither
# are likely not cycling and in G1 phase.
# 
# We assign scores in the CellCycleScoring() function, which stores S and G2/M scores in object meta data, 
# along with the predicted classification of each cell in either G2M, S or G1 phase. CellCycleScoring() 
# can also set the identity of the Seurat object to the cell-cycle phase by passing set.ident = TRUE 
# (the original identities are stored as old.ident). Please note that Seurat does not use the discrete
# classifications (G2M/G1/S) in downstream cell cycle regression. Instead, it uses the quantitative
# scores for G2M and S phase. However, we provide our predicted classifications in case they are of interest.

marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow[[]])

#Visualize the distribution of cell cycle markers across
RidgePlot(marrow, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)

# Regress out cell cycle scores during data scaling
# 
# We now attempt to subtract (‘regress out’) this source of heterogeneity from the data. 
# For users of Seurat v1.4, this was implemented in RegressOut. However, as the results of this 
# procedure are stored in the scaled data slot (therefore overwriting the output of ScaleData()),
# we now merge this functionality into the ScaleData() function itself.
# 
# For each gene, Seurat models the relationship between gene expression and the S and G2M cell 
# cycle scores. The scaled residuals of this model represent a ‘corrected’ expression matrix, that
# can be used downstream for dimensional reduction.

marrow <- ScaleData(marrow, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(marrow))
names(marrow)

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)

DimPlot(marrow)

# Alternate Workflow
# 
# The procedure above removes all signal associated with cell cycle. In some cases,
# we’ve found that this can negatively impact downstream analysis, particularly in
# differentiating processes (like murine hematopoiesis), where stem cells are quiescent 
# and differentiated cells are proliferating (or vice versa). In this case, regressing 
# out all cell cycle effects can blur the distinction between stem and progenitor cells as well.
# 
# As an alternative, we suggest regressing out the difference between the G2M and S phase
# scores. This means that signals separating non-cycling cells and cycling cells will be
# maintained, but differences in cell cycle phase among proliferating cells (which are 
# often uninteresting), will be regressed out of the data