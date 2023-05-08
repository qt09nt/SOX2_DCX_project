library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(ggplot2)
library(data.table)
library(tidyr)
library(Matrix)
library(ggrepel)
library(EnhancedVolcano)


if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

setwd("C:/Users/Diamandis Lab II/Documents/Queenie")

source("ramos_functions.R")

#install.packages("enrichR")

#March 28 2023 
#try implementing chooseR tool, which is a tool designed to work with Seurat for finding optimal number of clusters, for cluster stability
#gitub.com/rbpatt2019/chooseR
install.packages("renv")

renv::init(force = TRUE)


#remotes::install_github("wjawaid/enrichR")
library(enrichR); setEnrichrSite("Enrichr")

#load sample names for prenatal samples from Ramos dataset
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/late prenatal human neurodevelopment resolved by single nucleus transcriptomics/")
setwd("C:/Users/Diamandis Lab II/Documents/Queenie")

setwd("C:/Users/Diamandis Lab II/Documents/Queenie/ramos/GSM6720853")

ramos<- read.csv("samples_all.csv")
id<-ramos$id
prenatal<-id[1:30]
samples <-prenatal

#setwd("/data/gene_expression")

setwd("D:/")

project.name<-"ramos" 

dir.raw<-file.path(project.name)


sample.curr <- "GSM6720853"


# sample.curr<-samples[1]
# 
for(sample.curr in samples){

  setwd(paste0("D:/ramos/", sample.curr))
    # Load the PBMC dataset
    #pbmc.data <- Read10X(data.dir = "data/tmp/", gene.column=1)
    #pbmc.data <- Read10X(data.dir = "data/tmp/filtered_gene_bc_matrices/hg19", gene.column=1)
    #pbmc.data <- Read10X(data.dir = "C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/late prenatal human neurodevelopment resolved by single nucleus transcriptomics/filtered_gene_bc_matrices/hg19", gene.column=1)
    #pbmc.data <- Read10X(data.dir = getwd(), gene.column=1)
    pbmc.data <- Read10X(data.dir = getwd())
  
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = "ramos", min.cells = 3, min.features = 200)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    
    pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)
    
    
    # Visualize QC metrics as a violin plot
    #VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    
    
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
    
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)
    
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
    
 #optional
    DimPlot(pbmc, reduction = "pca", label=TRUE)
    ElbowPlot(pbmc)
    
    pbmc <- FindNeighbors(pbmc, dims = 1:12)
    #pbmc <- FindClusters(pbmc, resolution = 0.3)
    
  #optional for larger datasets increase the resolution value:  #changed from 2 to 4 for gsm6720853
    
    #(April 12 2023)
    #try out different values for resolution of gsm6730853 
    pbmc <- FindClusters(pbmc, resolution = 1.5)
    pbmc <- RunUMAP(pbmc, dims = 1:12)
   
      
  #optional 
    DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 5)
    VlnPlot(pbmc, features = "FABP7")
    
    #plot FABP7 expression over UMAP plot; 
    #adjust scale colours so none or low expression is closer to white;
    #high expression is blue
    FeaturePlot(pbmc, features = c("FABP7"), cols =c("darkgreen", "red"))
    
    
  #April 19 2023 save progress
    saveRDS(pbmc, "GSM6720853.rds")
    pbmc<-readRDS("GSM6720853.rds")
    
   
    ## find markers for every cluster compared to all remaining cells, 
    #report only the positive # ones
    pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
    #save the cluster marker genes as RDS file
    saveRDS(pbmc.markers, "GSM6720853_markers.rds")
    
    #heatmap for top 5 genes per cluster
    pbmc.markers %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = avg_log2FC) -> top5
    #DoHeatmap(pbmc, features = top5$gene) + NoLegend()
    
    #show legend, and change angle of text above the colour bar
    DoHeatmap(pbmc, features = top5$gene, angle = 45, size = 3) 
    
    #save pdf heatmap top 5 gene markers per cluster
    pdf(file = paste0("C:/Users/Diamandis Lab II/Documents/Queenie/ramos/output/March 23 2023/GSM6720853/", sample.curr, "_top5_gene_markers.pdf"),
        width = 12, height = 8)
    
############################# Cell Type Annotations    
    #The list of 500 most differentially expressed genes were calculated for each cluster from lines 331 to 
    #lines 338.
    #This list of genes for each cluster was then input into Enrichr at https://maayanlab.cloud/Enrichr/
    #Cell Type Annotations were retrieved from the Cell Types tab results page for the different databases
    
    #April 25, 2023 try labelling with one databases' results at a time
   #PanglaoDB Augmented 2021
    #April 27 2023 EDIT: PanglaoDB annotations are not as relevant for this dataset
    #since Purkinje Fiber Cells are part of the Cerebellum and this sample is taken from the
    #germinal matrix
    new.cluster.ids.panglaoDB<-c("Purkinje Fiber Cells",	
                                 "Purkinje Fiber Cells",	
                                 "Purkinje Fiber Cells",	
                                 "Purkinje Fiber Cells",	
                                 "Purkinje Fiber Cells",	
                                 "Immature Neurons",	
                                 "Purkinje Fiber Cells",	
                                 "Purkinje Fiber Cells",	
                                 "Immature Neurons",	
                                 "Purkinje Fiber Cells",	
                                 "Immature Neurons",	
                                 "Purkinje Fiber Cells",	
                                 "Purkinje Fiber Cells",	
                                 "Immature Neurons",	
                                 "Purkinje Fiber Cell",	
                                 "Radial Glia Cells",	
                                 "Purkinje Fiber Cells",	
                                 "Immature Neurons",	
                                 "Purkinje Fiber Cells",	
                                 "Immature Neurons",	
                                 "Purkinje Fiber Cells",	
                                 "Purkinje Fiber Cells",	
                                 "Oligodendrocyte Precursor cell:Brain",	
                                 "Purkinje Fiber Cells",	
                                 "Radial Glia Cells",	
                                 "Purkinje Fiber Cells",	
                                 "Purkinje Fiber Cells"	
    )
    
    names(new.cluster.ids.panglaoDB)<-levels(pbmc)
    pbmc<-RenameIdents(pbmc, new.cluster.ids.panglaoDB)
    DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
    
    
    
    #Azimuth Cell Types 2021
    new.cluster.ids.azimuth <-c("Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Layer 2-3 Glutamatergic Neuron, Intratelencephalon-Projecting CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Layer 6 Glutamatergic Neuron, Corticothalamic-Projecting CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679", 
                                "Sncg+ GABAergic Neuron",
                                "Non-neuronal Cell", 
                                "Glutamatergic Neuron CL0000679",
                                "GABAergic Neuron CL0000617",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Glutamatergic Neuron CL0000679",
                                "Non-neuronal Cell",
                                "Glutamatergic Neuron CL0000679",
                                "Sst+ GABAergic Neuron")
    
    names(new.cluster.ids.azimuth)<-levels(pbmc)
    pbmc<-RenameIdents(pbmc, new.cluster.ids.azimuth)
    DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
    
    ####Perform Differential Expression Analysis between Glutamatergic neuron clusters where FABP7 expression is 
    #relatively high (0, 1, 2, 3, 4, 6, 7, 22, 26) with other clusters that are the same cell type Glutamatergic neurons
    #that have low FABP7 expression (11, 8, 16, 12, 18, 9, 20, 25, 19, 23, 13)
    #https://satijalab.org/seurat/articles/de_vignette.html  
    
    # Perform default differential expression tests
    # 
    # The bulk of Seurat’s differential expression features can be accessed through the FindMarkers() function. As a default,
    # Seurat performs differential expression based on the non-parametric Wilcoxon rank sum test. This replaces the previous 
    # default test (‘bimod’). To test for differential expression between two specific groups of cells, specify the ident.1 and
    # ident.2 parameters.
    
    cluster11.markers <- FindMarkers(pbmc, ident.1 = 11, ident.2 = c(0, 1,2, 3, 4, 7, 22, 26), min.pct =0.25)
    
    #try to see if can compare the group of low FABP7 clusters to the higher FABP7 clusters Glutamatergic neurons
    clusterFABP7_low.markers <- FindMarkers(pbmc, ident.1 = c(11, 8, 16, 12, 18, 9, 20, 25, 19, 23, 13, 5), ident.2 = c(0, 1, 2, 3, 4, 6, 7, 22, 26), min.pct =0.25)
    
    #differential expression analysis of low FABP7 clusters of glutamatergic neurons, report only positive  ones (only.pos = TRUE)
    clusterFABP7_low_pos.markers <- FindMarkers(pbmc, only.pos = TRUE, ident.1 = c(11, 8, 16, 12, 18, 9, 20, 25, 19, 23, 13, 5), ident.2 = c(0, 1, 2, 3, 4, 6, 7, 22, 26), min.pct =0.25)
    
    ####### to do: do some pre-filtering for the differential expression analysis
    
    #save DE analysis results between FABP7 low Glutamatergic neurons vs the high FABP7 glutamatergic neurons
    #unfiltered DE results
    saveRDS(clusterFABP7_low.markers, "clusterFABP7_low.markers_unfiltered.rds")
    saveRDS(clusterFABP7_low_pos.markers, "clusterFABP7_low_pos.markers_unfiltered.rds")
    
    
    
    #view results 
    head(clusterFABP7_low.markers)
    
    # Prefilter features or cells to increase the speed of DE testing
    # 
    # To increase the speed of marker discovery, particularly for large datasets, Seurat allows for pre-filtering of 
    # features or cells. For example, features that are very infrequently detected in either group of cells, or features
    # that are expressed at similar average levels, are unlikely to be differentially expressed. Example use cases of 
    # the min.pct, logfc.threshold, min.diff.pct, and max.cells.per.ident parameters are demonstrated below.
    
    # The results data frame has the following columns :
    #   
    # p_val : p_val (unadjusted)
    # avg_log2FC : log fold-change of the average expression between the two groups. 
    # Positive values indicate that the feature is more highly expressed in the first group.
    # pct.1 : The percentage of cells where the feature is detected in the first group
    # pct.2 : The percentage of cells where the feature is detected in the second group
    # p_val_adj : Adjusted p-value, based on Bonferroni correction using all features in the dataset.
    
    
    #try to see if can compare the group of low FABP7 clusters to the higher FABP7 clusters Glutamatergic neurons
    #filter features that have less than a two-fold change between the average expression of low FABP7 glutamatergic
    #neurons and the high FABP7 glutamatergic neurons
    #log2(2) = 1
    clusterFABP7_low_log2.markers <- FindMarkers(pbmc, ident.1 = c(11, 8, 16, 12, 18, 9, 20, 25, 19, 23, 13, 5), 
                                            ident.2 = c(0, 1, 2, 3, 4, 6, 7, 22, 26), min.pct =0.25, logfc.threshold = log(2) )
    
    #try to see if can compare the group of low FABP7 clusters to the higher FABP7 clusters Glutamatergic neurons
    #Pre-filter features that are detected at <25% frequency in either low FABP7 Glutamatergic neurons or higher FABP7 Glutamatergic neurons
    clusterFABP7_low.markers <- FindMarkers(pbmc, ident.1 = c(11, 8, 16, 12, 18, 9, 20, 25, 19, 23, 13, 5), ident.2 = c(0, 1, 2, 3, 4, 6, 7, 22, 26), min.pct =0.25)
    
    #compare high FABP7 expression glutamatergic neurons clusters with clusters of low FABP7 expression  glutamatergic neurons
    clusterFABP7_high.markers <- FindMarkers(pbmc, ident.1 =c(0, 1, 2, 3, 4, 6, 7, 22, 26), ident.2 = c(11, 8, 16, 12, 18, 9, 20, 25, 19, 23, 13, 5), min.pct =0.25)
    
    saveRDS(clusterFABP7_high.markers, "clusterFABP7_high.markers.rds")
    
    clusterFABP7_high.markers<-readRDS("clusterFABP7_high.markers.rds")
    
    #results of diff expr wilcoxon rank sum test
    write.csv(clusterFABP7_high.markers, "clusterFABP7_high.markers.csv")
    
    write.csv(clusterFABP7_high_markers_t_test, "clusterFABP7_high_markers_t_test.csv")
    
    #### try pathway analysis 
    #https://jackbibby1.github.io/SCPA/index.html
    
    colnames(clusterFABP7_low_log2.markers)
    #"p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj"
    
    #make a column for the genes
    clusterFABP7_low_log2.markers$gene <-row.names(clusterFABP7_low_log2.markers)
    
    
    saveRDS(clusterFABP7_low_log2.markers, "clusterFABP7_low_log2.markers.RDS")
    readRDS("clusterFABP7_low_log2.markers.RDS")
    
    options(ggrepel.max.overlaps = Inf)
    
    #try to do the differential express test using t-test instead of the default Wilcoxon Rank Sum
    clusterFABP7_high_markers_t_test<-FindMarkers(pbmc, ident.1 = c(0, 1, 2, 3, 4, 6, 7, 22, 26), ident.2 = c(11, 8, 16, 12, 18, 9, 20, 25, 19, 23, 13, 5), test.use = "t")
    
    
    #plot the differentially expressed genes on a volcano plot
    #https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
    
    volcano_ggplot(top_table = clusterFABP7_low_log2.markers, p_value_cutoff = 0.1, subtitle="") 
    
    #######try plotting with Enhanced Volcano Plot package instead
    EnhancedVolcano(clusterFABP7_low_log2.markers,
                    lab = rownames(clusterFABP7_low_log2.markers),
                    subtitle = "Ramos et al (2022) GSM6720853", 
                    title = "low FABP7 expression Glutamatergic Neuron Clusters vs high FABP7 Glutamatergic Neuron Clusters ",
                    FCcutoff = 0.5,
                    x = 'avg_log2FC',
                    y = 'p_val_adj')
    
    EnhancedVolcano(clusterFABP7_low.markers,
                    lab = rownames(clusterFABP7_low.markers),
                    subtitle = "Ramos et al (2022) GSM6720853", 
                    title = "low FABP7 expression Glutamatergic Neuron Clusters vs high FABP7 Glutamatergic Neuron Clusters ",
                    FCcutoff = 0.5,
                    pCutoff = 0.05,
                    x = 'avg_log2FC',
                    y = 'p_val_adj')
    
    #plot the volcano plot of diff expr between high FABP7 glutamatergic neuron clusters vs low FABP7 glutamatergic neuron clusters
    EnhancedVolcano(clusterFABP7_high.markers,
                    lab = rownames(clusterFABP7_high.markers),
                    subtitle = "Ramos et al (2022) GSM6720853", 
                    title = "high FABP7 expression Glutamatergic Neuron Clusters vs low FABP7 Glutamatergic Neuron Clusters ",
                    FCcutoff = 0.5,
                    pCutoff = 0.05,
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    ylab = bquote(~-Log[10] ~ italic("adjusted p-value")),  #label the y axis label as -log10(adjusted p-value)
                    ylim = c(-5, 400))
    
    #
    #plot differentially expressed genes between high FABP7 expression glutamatergic neuron clusters
    #and low FABP7 expression glutamatergic neuron clusters using Seurat FindMarkers() Wilcoxon rank sum test
    #on y axis plot the original p values 
    enhanced_volcano(diff_expr_table = clusterFABP7_high.markers)
  
    
    #plot differentially expressed genes between high FABP7 expression glutamatergic neuron clusters
    #and low FABP7 expression glutamatergic neuron clusters using Seurat FindMarkers() t-test
    enhanced_volcano(diff_expr_table = clusterFABP7_high_markers_t_test)
    
    #try to do pathway analysis on the FABP7 low expression glutamatergic neuron clusters and the 
    #high FABP7 expression glutamatergic neuron clusters
    #https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html
    
    #Visualize the expression of the differentially expressed genes which are high in the high FABP7 glutamatergic neuron clusters on UMAP plots
    FeaturePlot(pbmc, features = c("STMN1", "RPL13", "MARCKSL1", "MT-CYB", "HSP90AA1", "PSIP1", "SOX4"), cols =c("darkgreen", "red"))
    
    #Visualize the expression of the differentially expressed genes which are high in the low FABP7 glutamatergic neuron clusters on UMAP plots
    FeaturePlot(pbmc, features = c("UBE4B", "RERE", "RGS7", "MSRA", "RBFOX1", "LRP1B", "ANK2", "PBX1", "SSBP2", "LRRC4C", "LSAMP"), cols =c("darkgreen", "red"))
    
    features.low <-c("UBE4B", "RERE", "RGS7", "MSRA", "RBFOX1", "LRP1B")
      
    ## try out different visualization methods May 5 2023
    RidgePlot(pbmc, features = features.low, ncol = 2)
    
    #################Pathway Analysis
    ################## get list of genes upregulated in the FABP7 high expression glutamatergic neurons group
    #sort by log 2 FC 
    
    websiteLive <- getOption("enrichR.live")
    if (websiteLive) {
      listEnrichrSites()
      setEnrichrSite("Enrichr") # Human genes   
    }
    
    #filter for the high FABP7 Glutamatergic Neuron group differentially expressed genes based on
    #log2 Fold Change greater than 0.5 
    high_FABP7_genes <- clusterFABP7_high.markers[clusterFABP7_high.markers$avg_log2FC >= 0.5, ]
    low_FABP7_genes <- clusterFABP7_high.markers[clusterFABP7_high.markers$avg_log2FC <= -0.5, ]
    
    #filter the high FABP7 Glutamatergic Neuron group differentially expressed genes further, 
    #based on adjusted pvalue < 0.05
    high_FABP7_genes <- high_FABP7_genes[high_FABP7_genes$p_val_adj <= 0.05, ]
    low_FABP7_genes<- low_FABP7_genes[low_FABP7_genes$p_val_adj <= 0.05, ] 
    
    high_FABP7_genes_list<-row.names(high_FABP7_genes)
    low_FABP7_genes_list<-row.names(low_FABP7_genes)
    
    
    #select favourite databases
    dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
    
    #query Enrichr with the genes upregulated in the FABP7 high expression glutamatergic neurons clusters
    if (websiteLive) {
      enriched.highFABP7 <- enrichr(high_FABP7_genes_list, dbs)
    }
    
    saveRDS(enriched.highFABP7, "enriched.highFABP7.rds")
    
    #query Enrichr for genes upregulated in the low FABP7 expression glutamatergic neurons
    if (websiteLive) {
      enriched.lowFABP7 <- enrichr(low_FABP7_genes_list, dbs)
    }
    
    setwd("C:/Users/Diamandis Lab II/Documents/Queenie/ramos/output/March 23 2023/GSM6720853/differential expression")
    
    saveRDS(enriched.lowFABP7, "enriched.lowFABP7.rds")
    
    #View results table
    gsm6720853_high_fabp7_GOBP_enrichr <-if (websiteLive) enriched.highFABP7[["GO_Biological_Process_2015"]]
    saveRDS(gsm6720853_high_fabp7_GOBP_enrichr, "gsm6720853_high_fabp7_GOBP_enrichr.RDS")
    write.csv(gsm6720853_high_fabp7_GOBP_enrichr, "gsm6720853_high_fabp7_GOBP_enrichr.csv")
    
    gsm6720853_low_fabp7_GOBP_enrichr <-if (websiteLive) enriched.lowFABP7[["GO_Biological_Process_2015"]]
    saveRDS(gsm6720853_low_fabp7_GOBP_enrichr, "gsm6720853_low_fabp7_GOBP_enrichr.RDS")
    write.csv(gsm6720853_low_fabp7_GOBP_enrichr, "gsm6720853_low_fabp7_GOBP_enrichr.csv")
    
    
    
    #Plot Enrichr GOBP output for high FABP7 expression Glutamatergic Neurons
    if (websiteLive){
      plotEnrich(enriched.highFABP7[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
    }
    
  
    #Plot Enrichr GOMF for high FABP7 expression Glutamatergic Neurons
    if (websiteLive){
      plotEnrich(enriched.highFABP7[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
    }
    
    
    #Plot Enrichr GOCC for high FABP7 expression Glutamatergic Neurons
    if (websiteLive){
      plotEnrich(enriched.highFABP7[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
    }
    
    
######### 
    
    #Plot Enrichr GOBP output for low FABP7 expression Glutamatergic Neurons
    if (websiteLive){
      plotEnrich(enriched.lowFABP7[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
    }
    
    
    #Plot Enrichr GOMF for low FABP7 expression Glutamatergic Neurons
    if (websiteLive){
      plotEnrich(enriched.lowFABP7[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
    }
    
    
    #Plot Enrichr GOCC for low FABP7 expression Glutamatergic Neurons
    if (websiteLive){
      plotEnrich(enriched.lowFABP7[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
    }
    
    
    seurat.clust<-pbmc@meta.data$seurat_clusters
    
    expr<-pbmc@assays$RNA@scale.data
    
    #expr<-pbmc@assays$RNA@counts
    
    clust.avg<-AverageExpression(pbmc)$RNA
    
    
    clust.names<-unique(pbmc@meta.data$seurat_clusters)
    
    list.genes<-lapply(clust.names, function(cl){
      #expr.curr<-expr[,seurat.clust==cl]
      #m<-rowMeans(expr.curr)
      #m.cutoff<-as.numeric(sort(m,decreasing=TRUE)[300])
      #genes.sel<-rownames(expr.curr)[m>=m.cutoff]
      #genes.sel<-rownames(expr.curr)[m>=0.4]
      
      #genes.sel<-rownames(clust.avg)[clust.avg[,cl]>2]
      #genes.sel<-c(genes.sel,rep("",5000-length(genes.sel)))
      
      avg.cutoff<-sort(clust.avg[,cl],decreasing=TRUE)[500]
      genes.sel<-rownames(clust.avg)[clust.avg[,cl]>=avg.cutoff]
      
      return(genes.sel)
    })
    
    # fix the error: "Error in app$vspace(new_style$`margin-top` %||% 0) : 
    #attempt to apply non-function"
    #for the line: 
    #df.genes<-bind_cols(list.genes)
    
    # DT1 = data.table(A=1:3,B=letters[1:3])
    # DT2 = data.table(A=4:5,B=letters[4:5])
    # l = list(DT1,DT2)
    # rbindlist(l)
    #March 14
    list.genes<-lapply(list.genes, list)
    
    #df.genes<-bind_cols(list.genes)
    
    df.genes <- rbindlist(list.genes, idcol = TRUE)
    
    df.genes$.id<-as.factor(df.genes$.id)
    df.genes$num <-c(rep(1:500, length(list.genes)))
    
    #re-format so dataframe df.genes goes from long to wide
    df_wide<-reshape(df.genes, direction = 'wide',     
                     idvar = '.id', 
                     timevar = 'num'
    )
    df.genes <-t(df_wide)
    df.genes <- df.genes[2:nrow(df.genes),]
    colnames(df.genes)<-clust.names
    df.genes <- as.data.frame(df.genes)
    write.table(df.genes, file=paste0("C:/Users/Diamandis Lab II/Documents/Queenie/ramos/output/02_genes/",sample.curr,".csv"),quote=FALSE,row.names=FALSE,sep=",")
    
    saveRDS(clust.avg, file=paste0("C:/Users/Diamandis Lab II/Documents/Queenie/ramos/output/02_cluster_expr/",sample.curr,".rds"))
    
    #  })
    #}
    
    #find the list of all available databases from Enrichr
    dbs <- listEnrichrDbs()
    
    #view and select favourite databases. then query enrich with the list of the genes.
    list.db<-list("CellMarker_Augmented_2021"=c("Oligodendrocyte Precursor cell:Brain",       # **
                                                "Oligodendrocyte:Brain",
                                                "Microglial cell:Embryonic Prefrontal Cortex"),# **
                  "PanglaoDB_Augmented_2021"="T Cells",
                  "Azimuth_Cell_Types_2021"=c("Oligodendrocyte CL0000128",                     # **
                                              "Microglia / Perivascular Macrophage CL0000881",
                                              "Astrocyte CL0000127",                           # **
                                              "Glutamatergic Neuron CL0000679",                # **
                                              "Proliferating Macrophage CL0000235")
                  
    )
    
    
    list.enriched<-list()
    df.genes <- as.data.frame(df.genes)
    list.enriched[colnames(df.genes)]<-list(NULL)
    #query the selected databases with the list of genes from df.genes
    for(clust in colnames(df.genes)){
      list.enriched[[clust]] <- enrichr(df.genes[,clust,drop=TRUE], names(list.db))    
    }
    
    saveRDS(list.enriched, file=paste0("C:/Users/Diamandis Lab II/Documents/Queenie/ramos/output/02_enrichr/", sample.curr,".rds"))
    
    cell.type<-NULL
    for(clust in colnames(df.genes)){
      
      enriched<-bind_rows(list.enriched[[clust]])
      enriched<-enriched[enriched$Term %in% do.call("c",list.db),]
      
      sig.total<-do.call("c", list.db)
      sig.missing<-sig.total[!sig.total %in% enriched$Term]
      
      if(length(sig.missing)>0){
        for(pos in 1:length(sig.missing)){
          enr.row<-tibble("Term"=sig.missing[pos], "Overlap"="", "P.value"=999, "Adjusted.P.value"=999, "Old.P.value"=999, "Old.Adjusted.P.value"=999,
                          "Odds.Ratio"=0, "Combined.Score"=0, "Genes"="")
          enriched<-rbind(enriched, enr.row)
        }
      }
      
      
      score.oligo_prec=enriched[enriched$Term=="Oligodendrocyte Precursor cell:Brain","Combined.Score"] # CellMarker_Augmented_2021
      
      score.oligo=enriched[enriched$Term=="Oligodendrocyte CL0000128","Combined.Score"] 
      
      score.astro=enriched[enriched$Term=="Astrocyte CL0000127","Combined.Score"]
      
      score.micro=enriched[enriched$Term=="Microglia / Perivascular Macrophage CL0000881","Combined.Score"]
      
      score.neuro=enriched[enriched$Term=="Glutamatergic Neuron CL0000679","Combined.Score"]
      
      
      score.mono_prolif=enriched[enriched$Term=="Proliferating Macrophage CL0000235","Combined.Score"]
      
      score.t_cell<-enriched[enriched$Term=="T Cells","Combined.Score"]
      
      # normally if oligo_prec > 300 the rest are quite low (<300)
      # but it can happen that oligo_prec=1700 and astro=1025; therefore, we also consider the spread
      is.gbm<-ifelse(score.oligo_prec>300 & 
                       score.oligo_prec>1.5*score.oligo & 
                       score.oligo_prec>1.5*score.neuro & 
                       score.oligo_prec>1.5*score.astro, TRUE, FALSE)
      
      is.oligo<-ifelse(score.oligo>300 & !is.gbm, TRUE,FALSE)
      
      is.astro<-ifelse(score.astro>300 & !is.gbm, TRUE,FALSE)
      
      is.micro<-ifelse(score.micro>300 & !is.gbm,TRUE,FALSE)
      
      is.neuro<-ifelse(score.neuro>300 & !is.gbm,TRUE,FALSE)
      
      is.endot<-FALSE
      
      
      
      is.mono_prolif<-ifelse(score.mono_prolif>300,TRUE,FALSE)
      
      is.t_cell<-ifelse(score.t_cell>100,TRUE,FALSE) # PanglaoDB signature, smaller values than Azimuth signatures
      
      
      row<-data.frame(oligodendrocyte=is.oligo,
                      astrocyte=is.astro,
                      microglia=is.micro,
                      neuron=is.neuro,
                      endothelial=is.endot,
                      glioblastoma=is.gbm,
                      monocyte_prolif=is.mono_prolif,
                      T_cell=is.t_cell)
      rownames(row)<-clust
      cell.type<-rbind(cell.type,row)
    }
    
    df.genes <-as.data.frame(df.genes)
    new.cluster.ids <-sapply(sort(names(df.genes)), function(clust){
      is.cell.type<-t(cell.type[rownames(cell.type)==clust,])[,1]
      if(sum(is.cell.type)==1){
        cell.name<-colnames(cell.type)[is.cell.type]
      }
      if(sum(is.cell.type)==0){
        cell.name<-"not_assigned"
      }
      if(sum(is.cell.type)>1){
        cell.name<-"multi_hit"
      }
      return(cell.name)
    })
    #names(new.cluster.ids) <- levels(pbmc)
    pbmc <- RenameIdents(pbmc, new.cluster.ids)
    
    p<-DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
    p<-p + labs(title=sample.curr)
    
    p <- p + scale_color_manual(breaks = c("oligodendrocyte", "astrocyte", 
                                           "microglia", "neuron", 
                                           "endothelial", "mono_prolif",
                                           "glioblastoma", "T_cell",
                                           "not_assigned", "multi_hit"),
                                values=c("royalblue", "turquoise", 
                                         "firebrick1","coral",
                                         "black","purple",
                                         "springgreen3","violet",
                                         "grey50","grey80"))
    print(p)
    
    png(filename=paste0("C:/Users/Diamandis Lab II/Documents/Queenie/ramos/output/02_plot/",sample.curr,".png"))
    print(p)
    dev.off()
    
    counts<-sapply(names(new.cluster.ids), function(clust.curr){
      sum(seurat.clust==clust.curr)
    })
    
    df.summary<-data.frame(cluster=names(new.cluster.ids),
                           cell_type=new.cluster.ids,
                           cell_count=counts)
    df.summary<-df.summary[order(as.numeric(df.summary$cluster)),]
    rownames(df.summary)<-NULL
    saveRDS(df.summary, file=paste0("C:/Users/Diamandis Lab II/Documents/Queenie/ramos/output/02_cell_type/",sample.curr,".rds"))
  #}) #try
  rm(list = ls())
    
  
}



