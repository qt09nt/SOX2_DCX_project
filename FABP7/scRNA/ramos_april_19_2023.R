library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(ggplot2)
library(data.table)
library(tidyr)
library(Matrix)

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
    FeaturePlot(pbmc, features = c("FABP7"), cols =c("white", "blue"))
    
    #optional do a heatmap for top 10 genes per cluster
    pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
    pbmc.markers %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = avg_log2FC) -> top5
    DoHeatmap(pbmc, features = top5$gene) + NoLegend()
    
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



