#boxplot for visualizing single gene for different time points
#data is grouped by time points (includes all 3 cell lines)

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/prot_org_secr")

library(viridis)

annotations <- read.csv("Annotations.csv")
working_data <- read.csv("Working_data.csv")

row.names(working_data) <- working_data$genes
working_data[,1]<-NULL

working_data$genes <- NULL

#get rid of the first column
annotations[,1] <- NULL

#For annotations, set sample name as row.names
row.names(annotations)<-annotations$sample 
annotations$timepoints <- as.factor(annotations$timepoints)


library(ggplot2)


get_gene_expr<-function(gene, working_data){
  gene_expression <- working_data[gene, ]
  gene_expr_t <- t(gene_expression)
  #change gene_expr_t row names to match the row names of annotation file
  row.names(gene_expr_t) <-row.names(annotations)
  #merge the single gene expression data with the annotations file 
  combined <-cbind(annotations, gene_expr_t)
  
  #re-order timepoints
  
  combined$timepoints <<- ordered(combined$timepoints, levels = c("Week 4", "Week 6", "Week 8", "Week 10",
                                                                "Week 12", "Week 15", "Week 16"))
}

#gene is gene name of interest
#combined is the annotations dataframe combined with gene expression data for 
#that gene


ggplot(combined, aes(x=timepoints, y=CLCA1)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(paste0(" Gene Expression Over Time")) +
    xlab("timepoints")
    ylab("Gene expression")




