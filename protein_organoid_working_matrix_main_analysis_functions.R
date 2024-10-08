# Helper functions for main_analysis

get_gene_names = function(data, gene_colname = "Gene names", protein_colname = "Protein names"){
  # go through the gene names and pick the first one if there are multiple
  # if there are no gene names, take the protein name
  # expects data to be perseus output with columns 'Gene names' and 'Protein names'
  
  multi_names = strsplit(data[[gene_colname]], ";")
  protein_multi_names = strsplit(data[[protein_colname]], ";")
  gene_names = rep(NA, nrow(data))
  for(i in 1:length(multi_names)){
    gene_names[i] = multi_names[[i]][1]
    if(is.na(gene_names[i])){
      gene_names[i] = protein_multi_names[[i]][1]
    }
  }
  return(make.names(gene_names, unique = T))
}

# gets the top differentially expressed proteins using the limma package
#see limma reference manual p.233 for topTable function
#logFC = log2 fold change of exon vs other exons for the same gene
#t = moderated t-statistic
#F = moderated F statistic
#P.Value = p-value
#FDR = false discovery rate 
#looks like p-values are adjusted by fdr
#

get_top_table = function(data, design){
  fit = lmFit(data, design) 
  fit = eBayes(fit) 
  top_table =  topTable(fit, coef="difference", adjust="fdr", number = nrow(data))
  return(top_table)
}

# helper for getting differentially expressed proteins at specific timepoint
get_tp_top_table = function(data, groups, cell_type, timepoint){
  cell_type_data = data[, groups$cell_type == cell_type]
  tps = groups$tp[groups$cell_type == cell_type]
  
  difference = rep(0, ncol(cell_type_data))
  difference[tps == timepoint] = 1
  design = cbind(control=1, difference=difference)
  
  top_table =  get_top_table(cell_type_data, design)
  return(top_table)
}

# #function for getting differentially expressed proteins for different cultures with specific media
# get_culture_top_table = function(data, groups, culture, media){
#   culture_data = data[, groups$culture == culture]
#   medias = groups$media_num[groups$culture == culture]
# 
#   difference = rep(0, ncol(culture_data))
#   difference[medias == media] = 1
#   design = cbind(control=1, difference = difference)
# 
#   top_table = get_top_table(culture_data, design)
#   return(top_table)
# }


#volcano plot where p. value cut off is 0.1
volcano_plot = function(top_table, title){
  v = EnhancedVolcano(top_table, x="logFC", y= "adj.P.Val", lab=row.names(top_table), title = title,
                      pCutoff = 0.1, FCcutoff = 1.333, 
                      ylim = c(0, max(-log10(top_table$adj.P.Val), na.rm=TRUE)),
                      col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
                      colAlpha = 1.0)
  #pdf(paste("../figures/", title, "volcano.pdf"), width = 8, height = 8)
  pdf(paste("figures/timepoints comparison/", title, "volcano.pdf"), width = 8, height = 8)
  print(v)
  dev.off()
  return(v)
}

#function for comparing timepoints while including all cell types
#annotations is the meta data file, working_data is the working matrix processed from Perseus
#timepoint1 is the first time point to compare; timepoint2 is the second timepoint to compare
get_timepoints <- function(annotations, working_data, timepoint1, timepoint2){
  tp_meta <- annotations[annotations$timepoint == timepoint1 | annotations$timepoints == timepoint2, ]
  
  #get specific samples for the 2 time points
  subset_samples <-c(tp_meta$sample)          #copy sample names into a vector
  subset_expr <- working_data %>%             #subset just the sample columns
    select(subset_samples)
  difference = ifelse(tp_meta$timepoint == timepoint1, 1, 0)
  design = cbind(control=1, difference = difference)
  tp1_vs_tp_2_top_table = get_top_table(subset_expr, design)
  
  #save results
  write.table(tp1_vs_tp_2_top_table, paste0("results/", timepoint1, "_vs_", timepoint2, "_top_table.txt"), sep = "\t", quote = F, row.names = T, col.names = T)
  print(volcano_plot(tp1_vs_tp_2_top_table, paste0("Organoid Protein Secretome - All Cell Types - ", timepoint1, " vs ", timepoint2) ))
  return(tp1_vs_tp_2_top_table)
}

#function for individual cell types time point comparisons
get_timepoints_individual_cell_types <- function(annotations, working_data, timepoint1, timepoint2, cell_type){
  tp_meta <- annotations[annotations$timepoint == timepoint1 | annotations$timepoints == timepoint2, ]
  
  #get specific samples for the 2 time points
  subset_samples <-c(tp_meta$sample)          #copy sample names into a vector
  subset_expr <- working_data %>%             #subset just the sample columns
    select(subset_samples)
  difference = ifelse(tp_meta$timepoint == timepoint1, 1, 0)
  design = cbind(control=1, difference = difference)
  tp1_vs_tp_2_top_table = get_top_table(subset_expr, design)
  
  #commented out the saving top table part and volcano plot because already have those saved; uncomment when needed
  # write.table(tp1_vs_tp_2_top_table, paste0("results/SYN_only/", cell_type, "_", timepoint1, "_vs_", timepoint2, "_top_table.txt"), sep = "\t", quote = F, row.names = T, col.names = T)
  # print(volcano_plot(tp1_vs_tp_2_top_table, paste0("Organoid Protein Secretome -", cell_type, "_", timepoint1, " vs ", timepoint2) ))
  return(tp1_vs_tp_2_top_table)
}

#volcano plot where p value cutoff is 0.05
volcano_plot = function(top_table, title){
  v = EnhancedVolcano(top_table, x="logFC", y= "adj.P.Val", lab=row.names(top_table), title = title,
                      pCutoff = 0.05, FCcutoff = 1.333, 
                      ylim = c(0, max(-log10(top_table$adj.P.Val), na.rm=TRUE)),
                      col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
                      colAlpha = 1.0)
  pdf(paste("../figures/", title, "volcano.pdf"), width = 8, height = 8)
  print(v)
  dev.off()
  return(v)
}


#function for taking a dataframe of GSEA enriched annotations (ie. gsea_report_for_na_pos_1641334710091.tsv) 
#function requires as input the GSEA pre-ranked results (both positive NES (enriched in the first timepoint) and negative
#NES (enriched in the samples in the second timepoint being compared)
#which contains the enriched pathways and creates a plot of the enrichment scores between 2 different timepoints
#output is the scatterplot of the pathway annotations and their FDR q. value and NES score

library("ggplot2")
library("dplyr")
library("viridis")

# #Plot FDR.q.val (different colour scale), and NES score (different size of dots)
#  <- function(gsea_results)

#function to process the GSEA pre-ranked results & plot as scatterplot
#input gsea_report_na_pos is the filename for the positive gsea report results (showing enrichment of pathways in first timepoint )
#timepoint1 is the string name of timepoint 1 ie. "week 4"

#for fdr q value cut off, adjust as necessary:
# If there aren't any interesting pathways
# you should use less stringent cut off
# if there are a lot of pathways you should use more stringent
# the idea is to make a good story
# remember that statistical significance is not necessarily representative of biological importance, take a transcription factor a small change in abundance can have very small significance statistically but a very big change biologically due to amplification steps in transcription, translation and post translational modificiation

#https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_False_Discovery_Rate

#function for filtering GSEA results for DCX and SOX2 
filter_gsea <- function(gsea_report_na_pos, gsea_report_for_na_neg, cell_type1, cell_type2){
  
  #fill = TRUE parameter is added to add "NA" values if there are missing cells within the GSEA file result
  DCX_gsea <-read.table(file = gsea_report_na_pos, sep = '\t', header = TRUE, fill = TRUE)
  colnames(DCX_gsea)[1] <-"Pathways"     #rename first column of GSEA pos results to "Pathways"
  #filter the GSEA results dataframe to keep just the pathways which have a NOM p-value of less than 0.07 and a 
  #FDR q values less than 0.1
  filtered_DCX <- DCX_gsea[DCX_gsea$`NOM.p.val` <= 0.05, ]
  filtered_DCX <- filtered_DCX[filtered_DCX$`FDR.q.val` <= 0.01, ]
  filtered_DCX <- filtered_DCX[,1:11]
  filtered_DCX$cell_type <- cell_type1
  filtered_DCX$NES <- as.numeric(filtered_DCX$ES)
  filtered_DCX$NOM.p.val<- as.numeric( filtered_DCX$NOM.p.val)
  filtered_DCX <<- filtered_DCX
  
  # #load in the gsea report for na_neg (negative enrichment for DCX compared to SOX2)
  SOX2_gsea <- read.table(file = gsea_report_for_na_neg, sep = '\t', header = TRUE, fill=TRUE)
  colnames(SOX2_gsea)[1] <- "Pathways"   ##rename first column of GSEA pos results to "Pathways"
  filtered_SOX2 <- SOX2_gsea[SOX2_gsea$`NOM.p.val` <= 0.1,]
  filtered_SOX2 <- filtered_SOX2[filtered_SOX2$`FDR.q.val` <= 0.1, ]
  filtered_SOX2 <- filtered_SOX2[,1:11] #keep columns 1 to 12 only
  filtered_SOX2$cell_type <- cell_type2
  filtered_SOX2$NES <- as.numeric(filtered_SOX2$ES)
  filtered_SOX2$NOM.p.val<- as.numeric( filtered_SOX2$NOM.p.val)
  filtered_SOX2 <<- filtered_SOX2 
}  
  
# function to plot the NES scores for DCX GSEA and SOX2 GSEA pathways on separate plots
#cell_type is a string specifying cell type ie. "DCX" or "SOX2"
#filtered_df is the filtered gsea results dataframe (filtered by NOM pvalue and FDR q value)
plot_NES <-function(cell_type, filtered_df){
  p <- ggplot(filtered_df, aes(x = cell_type, y = Pathways, size = NES, color = FDR.q.val )) +
    geom_point(alpha = 0.7) + #scatterplot
    scale_size(range = c(2, 5), name = "NES") +  #change the size of the points and create a name for the size legend
    scale_color_viridis(discrete = FALSE)
  
  plot(p)
  #function for saving figures as pdf, png or jpg
  #set directory for saving GSEA results comparing organoid secretome timepoints
  #setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/figures/timepoint comparisons GSEA")
  
  #set directory for GSEA results for comparing organoid DCX vs SOX2 pathways
  setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/figures/")
  
  
  #save as png file, to save as pdf write "pdf", or jpg, write "jpg"; res = resolution
  png(filename=paste0(cell_type, "_gsea_plot.pdf"), width=3000, height=3000, res=300)
  
  #function that makes the plot
  plot(p)
  dev.off()
}  

#function for plotting merged SOX2 & DCX gsea pathways 
plot_NES_merged <- function(filtered_df){
  
   p <- filtered_df %>%
    mutate(Pathways = fct_reorder(Pathways, desc(FDR.q.val))) %>%
    ggplot(aes(x = cell_type, y = Pathways, size = NES, color = FDR.q.val )) +
      geom_point(alpha = 0.7) + #scatterplot
      scale_size(range = c(2, 5), name = "NES") +  #change the size of the points and create a name for the size legend
      scale_color_viridis(discrete = FALSE)
  
  plot(p)
  #function for saving figures as pdf, png or jpg
  #set directory for saving GSEA results comparing organoid secretome timepoints
  #setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr/figures/timepoint comparisons GSEA")
  
  #set directory for GSEA results for comparing organoid DCX vs SOX2 pathways
  setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/figures/")
  
  
  #save as png file, to save as pdf write "pdf", or jpg, write "jpg"; res = resolution
  #png(filename=paste0( "DCX_vs_SOX2_GOBP_gsea_plot.pdf"), width=3000, height=3000, res=300)
  
  #width and height input is in inches
  pdf(file="DCX_vs_SOX2_GOBP_gsea_plot.pdf", width=13, height=13)
  
  
  #function that makes the plot
  plot(p)
  dev.off()
}  

  
#function written by Brian Lam
#used to get the NOM.P.val < cutoff value
process.csv <- function(x,filter=TRUE,cutoff=0.1){
  #'process.csv
  #'
  #'This function takes a csv file with separators values and creates a dataframe
  #'
  #'@param x this is the input for a csv file to be converted into a dataframe
  #'@param cutoff this is an optional cutoff value. The default is 0.1.
  #process the csv file
  csvfile <- read.table(x,header=TRUE, sep = '\t', fill = TRUE)
  #create a vector using booleans to identify the elements of NOM.p.val that are < cutoff.
  if (filter == TRUE) {
    boolean.greater <- csvfile[(csvfile$`NOM p-val`<cutoff),]
    return(boolean.greater)
  } else {
    return(csvfile)
  }
  #returns the complete dataframe with the values spliced out according to the cutoff.
}


process.csv2 <- function(x,filter=TRUE,FDRqvalcutoff=0.15){
  #'process.csv
  #'
  #'This function takes a csv file with separators values and creates a dataframe
  #'
  #'@param x this is the input for a csv file to be converted into a dataframe
  #'@param cutoff this is an optional cutoff value. The default is 0.1.
  #process the csv file
  csvfile <- read.table(x,header=TRUE, sep = '\t', fill = TRUE)
  #keep just the pathways (rows) where NOM p. value < 0.1, and FDR q.value < 0.15
  filtered <- csvfile[csvfile$`NOM.p.val` <= 0.07, ]
  filtered <- filtered[filtered$`FDR.q.val` <= FDRqvalcutoff, ]
  #return(filtered)
}

#function to save processed/filtered GSEA files:
save_processed_csv <- function(gsea_pos, gsea_neg, timepoint1, timepoint2, processed_results_directory, timepoint_comparisons){
  tp1 <- as.data.frame(process.csv2(gsea_pos))
  tp2 <- as.data.frame(process.csv2(gsea_neg))
  results.folder <-dir.create(paste0(processed_results_directory, timepoint_comparisons))
  setwd(paste0(processed_results_directory, timepoint_comparisons))
  write.csv(tp1, paste0(timepoint1, ".csv"))
  write.csv(tp2, paste0(timepoint2, ".csv"))
            
}

process.csv <- function(x,filter=TRUE,cutoff=0.1){
  #'process.csv
  #'
  #'This function takes a csv file with separators values and creates a dataframe
  #'
  #'@param x this is the input for a csv file to be converted into a dataframe
  #'@param cutoff this is an optional cutoff value. The default is 0.1.
  #process the csv file
  csvfile <- read.csv(x,header=TRUE)
  #create a vector using booleans to identify the elements of NOM.p.val that are < cutoff.
  if (filter == TRUE) {
    boolean.greater <- csvfile[(csvfile[,"NOM.p.val"]<cutoff),]
    return(boolean.greater)
  } else {
    return(csvfile)
  }
  #returns the complete dataframe with the values spliced out according to the cutoff.
}

unique.comparisons <- function(dataframe1,dataframe2){
  #'unique.comparisons
  #'
  #'This function takes two vectors as inputs and returns the elements in
  #'vector1 that are unique to both cases
  #'
  #'@param dataframe1 this is the first vector to be compared to the second vector1
  #'@param dataframe2 this is the second vector to be compared
  #splice out the names of the processes or whatever is being compared.
  #this was originally designed to compare targets in GSEA
  vector1 <- dataframe1[,'Pathways']
  vector2 <- dataframe2[,'Pathways']
  #generates a vector of booleans that show elements that are unique to vector1
  #when compared to each of the other vectors
  # boolean.comparison1 <- (!(vector1 %in% vector2) & !(vector1 %in% vector3)
  #                         & !(vector1 %in% vector4) & !(vector1 %in% vector5))
  
  boolean.comparison1 <- (!(vector1 %in% vector2))
  #gets only the rows that are in unique from the original dataframe of vector1
  unique.dataframe <- dataframe1[boolean.comparison1,]
  #returns the dataframe for future use
  return(unique.dataframe)
}

common.comparisons <- function(dataframe1,dataframe2){
  #'Common.comparisons
  #'
  #'This function takes two vectors as inputs and returns the elements in
  #'vector1 that are common to both cases
  #'
  #'@param dataframe1 this is the first vector to be compared to the second vector1
  #'@param dataframe2 this is the second vector to be compared
  #splice out the names of the processes or whatever is being compared.
  #this was originally designed to compare targets in GSEA
  vector1 <- dataframe1[,'Pathways']
  vector2 <- dataframe2[,'Pathways']
  #gets only the rows that are in common
  common.dataframe <- dataframe1[vector1 %in% vector2,]
  #returns the dataframe for future use
  return(common.dataframe)
}

#function to fill out diffexpressed top table, which is used for labelling volcano plots
#pvalue can be 0.05 or 0.10
get_diffexpressed <- function(top_table, p_value, condition1, condition2){
  #if log2FoldChange > 1.5 and pvalue <0.1 set as "UP in GLICO"
  top_table$diffexpressed[top_table$logFC > 1.5 & top_table$adj.P.Val < p_value] <- paste0("UP in ", condition1)
  
  #if log2FoldChange < -1.5 and pvalue < 0.1, set as "UP in monolayer"
  top_table$diffexpressed[top_table$logFC < -1.5 & top_table$adj.P.Val < p_value] <- paste0("UP in ", condition2)
  
  #label extra gene candidates which have a log2FoldChange > 3 or < -3 for potential targets in downstream experiments
  top_table$diffexpressed[top_table$logFC < - 3 & top_table$adj.P.Val > p_value ] <- "< -3 log2Fold Change"
  top_table$diffexpressed[top_table$logFC > 3 & top_table$adj.P.Val > p_value ] <- "> 3 log2Fold Change"
  
  # Now write down the name of genes beside the points...
  # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  top_table$delabel <- NA
  top_table$delabel[top_table$diffexpressed != "NO"] <- top_table$gene[top_table$diffexpressed != "NO"]
  
  top_table <<- top_table
}

#function to fill out diffexpressed top table, but don't highlight <-3 log2FC and >3 log2FC:
get_diffexpressed_modified <- function(top_table, p_value, condition1, condition2){
  #if log2FoldChange > 1.5 and pvalue <0.1 set as "UP in GLICO"
  top_table$diffexpressed[top_table$logFC > 1.5 & top_table$adj.P.Val < p_value] <- paste0("UP in ", condition1)
  
  #if log2FoldChange < -1.5 and pvalue < 0.1, set as "UP in monolayer"
  top_table$diffexpressed[top_table$logFC < -1.5 & top_table$adj.P.Val < p_value] <- paste0("UP in ", condition2)
  
  #if log2FoldChange < -1.5 and pvalue > 0.1 set up as Log2FC:
  top_table$diffexpressed[top_table$logFC < -1.5 & top_table$adj.P.Val > p_value] <- paste0("Log2FC")
  top_table$diffexpressed[top_table$logFC > 1.5 & top_table$adj.P.Val > p_value] <- paste0("Log2FC")
  
  # Now write down the name of genes beside the points...
  # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  top_table$delabel <- NA
  top_table$delabel[top_table$diffexpressed != "NO"] <- top_table$gene[top_table$diffexpressed != "NO"]
  
  top_table <<- top_table
}




#function for plotting volcano plots using ggplot 
#p_value_cutoff can be 0.1, 0.05
volcano_ggplot <- function(top_table, p_value_cutoff, subtitle){
  ggplot(data=top_table, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +    
    scale_color_manual(values=c("#F68D91",  "#4E80A5", "grey", "#A4D49C", "orange")) +
    geom_vline(xintercept=c(-1.5, 1.5), col="red") +
    geom_hline(yintercept=-log10(p_value_cutoff), col="red") +
    labs(title = 'GLICOs dataset',
         subtitle=paste0(subtitle, " p. adj. value cut-off = ", p_value_cutoff),
         x = 'log2 Fold Change')
}

#version 2 - colour in data which pass log2FC cutoff
modified_volcano_ggplot <- function(top_table, p_value_cutoff, subtitle){
  ggplot(data=top_table, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel(    #if top_table$genelabels = TRUE, label the gene with the gene name from top_table$gene
      aes(logFC, -log10(adj.P.Val)),
      size = 5,
      label = ifelse(top_table$genelabels, top_table$gene, ""),
      box.padding = unit(0.45, "lines"),
      colour = "black",
      hjust = 1
    ) +    
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 18),  #legend font size
          plot.title = element_text(color = "black", size = 15, face = "bold"),
          plot.subtitle = element_text(color = "black", size = 12),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")) + 
    scale_color_manual(values=c("#A4D49C", "grey", "red", "blue","orange")) +
    geom_vline(xintercept=c(-1.5, 1.5), col="red") +
    labs(title = 'GLICOs dataset', 
         subtitle=paste0(subtitle, " p. adj. value cut-off = ", p_value_cutoff),
         x = 'log2 Fold Change')
}

#version 1 - keep data that pass log2FC cutoff as grey
modified_volcano_ggplot <- function(top_table, p_value_cutoff, subtitle){
  ggplot(data=top_table, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel(    #if top_table$genelabels = TRUE, label the gene with the gene name from top_table$gene
      aes(logFC, -log10(adj.P.Val)),
      size = 5,
      label = ifelse(top_table$genelabels, top_table$gene, ""),
      box.padding = unit(0.45, "lines"),
      colour = "black",
      hjust = 1
    ) +    
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 18),  #legend font size
          plot.title = element_text(color = "black", size = 15, face = "bold"),
          plot.subtitle = element_text(color = "black", size = 12),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")) + 
    scale_color_manual(values=c("grey", "red", "blue", "#A4D49C","orange")) +
    geom_vline(xintercept=c(-1.5, 1.5), col="red") +
    labs(title = 'GLICOs dataset', 
         subtitle=paste0(subtitle, " p. adj. value cut-off = ", p_value_cutoff),
         x = 'log2 Fold Change')
}

get_pscore_table = function(top_table){
  direction_mult = sign(top_table$logFC)
  signed_pscore = -log10(top_table$adj.P.Val) * direction_mult
  pscore_table = cbind.data.frame(row.names(top_table), signed_pscore)
  pscore_table = pscore_table[order(pscore_table$signed_pscore, decreasing = T),]
  return(pscore_table)
}

get_boxplot = function(data, groups, title, adj_mult = 1,
                       order_groups = c(), comparisons = list(),
                       show_pair_pval = F, font_size = 15, line_thickness = 1.0,
                       x_lab = NULL, y_lab = "Expression", jitter = T){
  
  df = cbind.data.frame(data = as.numeric(data), groups = groups)
  
  if(length(order_groups) > 0){
    df$groups = factor(df$groups, levels = order_groups)
  }
  
  if(jitter){
    add = "jitter"
  }else{
    add = NULL
  }
  
  g =  ggboxplot(df, x = "groups", y = "data", color = "groups", 
                 add = add, legend = "none", size = line_thickness,
                 palette = c("#9DB4DD", "#F69679", "#F7AFC3",  "#A4D49C", "#FDC98F",
                             "#AB92C5", "#4E80A5", "#BD7691", "#FFF79F")) +
    rotate_x_text(angle = 45) 
  
  if(!is.null(x_lab)){
    g = g + xlab(x_lab)
  }
  if(!is.null(y_lab)){
    g = g + ylab(y_lab)
  }
  
  if(length(unique(df$groups)) == 2){
    df$groups = as.character(df$groups)
    groups1 = df$data[df$groups == unique(df$groups)[1]]
    groups2 = df$data[df$groups == unique(df$groups)[2]]
    test = t.test(groups1, groups2)
    
    if(adj_mult == 1){
      tab = cbind.data.frame(test$p.value)
      colnames(tab) = c("p.value")
      p_tab <- tableGrob(tab, rows = NULL)
    }else{
      tab = rbind.data.frame(test$p.value, (test$p.value)*adj_mult)
      rownames(tab) = c("p", "adj.p")
      p_tab <- tableGrob(tab, cols = NULL)
    }
    
    g = g + geom_signif(comparisons = list(unique(df$groups)), 
                        test = "t.test",
                        map_signif_level=c("***"=(0.001/adj_mult), "**"=(0.01/adj_mult), "*"=(0.05/adj_mult)), 
                        step_increase = c(0.1, 0.1, 0.1),
    )
    
    if(show_pair_pval){
      g = ggdraw() + draw_plot(g, width = 0.7) + draw_plot(p_tab, x = 0.7, width = 0.3)
    }
  }else if(length(comparisons > 0)){
    g = g + geom_signif(comparisons = comparisons, 
                        test = "t.test",
                        map_signif_level=c("***"=(0.001/adj_mult), "**"=(0.01/adj_mult), "*"=(0.05/adj_mult)), 
                        step_increase = c(0.1, 0.1, 0.1), +
                          geom_hline(yintercept = mean(df$data), linetype = 2) # Add horizontal line at base mean
                        
    )
  }else{
    if(show_pair_pval){
      label_type = "p.format"
    }else{
      label_type = "p.signif"
    }
    
    g = g + stat_compare_means(method = "anova", label.y = (max(df$data) + sd(df$data)),
                               label.x.npc = 0.1)+ # Add global annova p-value
      stat_compare_means(label = label_type, method = "t.test",  # alternatively: label = "p.format"
                         ref.group = ".all.", hide.ns = T, label.y.npc = "top")+# Pairwise comparison against all
      geom_hline(yintercept = mean(df$data), linetype = 2) # Add horizontal line at base mean
  }
  
  # apply font size
  g = g + theme(text = element_text(size = font_size))
  g = g + ggtitle(title)
  
  return(g)
}

venn_diagram = function(list_of_sets, names, title, output=T){
  myCol <- brewer.pal(max(3, length(names)), "Pastel2")[1:length(names)]
  
  v = venn.diagram(
    x = list_of_sets,
    category.names = names,
    filename = paste("figures/", title, ".tiff"),
    output=output,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans"
  )
  return(v)
}




