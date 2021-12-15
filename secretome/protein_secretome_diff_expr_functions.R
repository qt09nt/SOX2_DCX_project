#Functions for analysis of Sofia's protein organoid secretome differential expression

#article on Why, When and How to Adjust Your P Values?
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099145/
#doi: 10.22074/cellj.2019.5992
#Jafari M, Ansari-Pour N. Why, When and how to adjust your P values? Cell J. 2019; 20(4): 604-607. doi: 10.22074/cellj.2019.5992.

#Adjust P-values for Multiple Comparisons
#https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html

#get columns for specific clusters and rename column names with cluster assignment
subset_cluster<- function(tsne_clusters, cluster.name){
  keep_sel <- tsne_clusters[,"cluster"] == cluster.name
  sub_cluster <- prot_organoid_secretome_log[,keep_sel]        #expression data for the immature cluster only 
  colnames(sub_cluster) <- paste(cluster.name, colnames(sub_cluster), sep = "_") #rename samples with cluster assignment
  cluster_df <- data.frame(sub_cluster)
  return(cluster_df)
}

#ANOVA testing differences between groups
doANOVA<-function(group,data.row){
  if(max(data.row, na.rm=TRUE) - min(data.row, na.rm=TRUE) != 0){
    names(group)<-names(data.row)
    data.sub<-data.frame(variable=data.row,group=group)
    colnames(data.sub)<-c("variable","group")
    data.sub<-data.sub[is.finite(data.sub$variable),]
    p.value<-summary(aov(variable~group,data=data.sub))[[1]]['group','Pr(>F)']
  }else{
    p.value<-999
  }
  return(p.value)
}

calcDiffExpr<-function(label.g1, label.g2){
  message("calcDiffExpr: ",label.g1,"|",label.g2)
  g1<-prot[,grepl(label.g1,colnames(prot))]
  g2<-prot[,grepl(label.g2,colnames(prot))]
  
  prot.sel<-cbind(g1, g2)
  group<-c(rep(label.g1, ncol(g1)), rep(label.g2, ncol(g2)))
  
  res_anova<-apply(prot.sel, 1, function(x){
    pval=doANOVA(group=group, data.row=x)
    return(data.frame(pval=pval))
  })
  
  res_logfc <- apply(prot.sel, 1, function(x){
    mean.g1=mean(x[group==unique(group)[1]])
    mean.g2=mean(x[group==unique(group)[2]])
    fold.change<-(2^mean.g1-2^mean.g2) / 2^mean.g1 # expression values in log scale
    return(data.frame(fold.change=fold.change))
  })
  
  
  df_anova <- data.frame(gene=rownames(prot.sel), bind_rows(res_anova))
  df_logfc <- data.frame(gene=rownames(prot.sel), bind_rows(res_logfc))
  # 
  df <- cbind(df_anova, df_logfc)
  
  #get rid of duplicate gene column
  df[,3]<-NULL
  
  #Bonferroni multiple testing adjustment
  #multiplies the raw P values by the number of tests m (i.e. length of the vector P_values). Using the p.adjust 
  #function and the 'method' argument set to "bonferroni", we get a vector of same length but with adjusted P values.
  #This adjustment approach corrects according to the family-wise error rate of at least one false positive
  #(FamilywiseErrorRate (FWER)=Probability (FalsePositive???1))
  
  df$pval.adj.bonferroni <- p.adjust(df$pval, method="bonferroni")
  
  #select for differentially expressed genes 
  df.sel <- df[df$pval<0.05 & abs(df$fold.change)>0.5,]
  
  
  #get rid of rows where the gene column is NA
  df.sel_no_NAs <- df.sel %>% drop_na(gene)
  
  ####adjust p.values to correct for multiple testing
  
  #Benjamini-Hochberg/FDR
  #false discovery rate (FalseDiscoveryRate (FDR)=Expected (FalsePositive/ (FalsePositive+TruePositive))). 
  #In other words, FDR is the expected proportion of false positives among all positives which rejected 
  #the null hypothesis and not among all the tests undertaken. In the FDR method, P values are ranked in an 
  #ascending array and multiplied by m/k where k is the position of a P value in the sorted vector and m is 
  #the number of independent tests. 
  
  df.sel_no_NAs$pval.adj.fdr <- p.adjust(df.sel_no_NAs$pval, method="fdr")
  df.sel_no_NAs <- df.sel_no_NAs
  
  #fname.out<-file.path(dir.results, "04_diff_gene_expression",paste0("diff_",label.g1,"_",label.g2,".csv"))
  
  #remove the genes which have "-Inf" in the fold.change column
  df_final <-df.sel_no_NAs[!(df.sel_no_NAs$fold.change == "-Inf"),]
  
  saveRDS(df_final, file = paste0("diff_expr_padj_",label.g1,"_",label.g2,".rds"))
  fname.out<-file.path(dir.results, paste0("diff_expr_padj_",label.g1,"_",label.g2,".csv"))
  
  write.csv(df_final, file=fname.out,row.names=FALSE, quote = FALSE)
}

#alternative calculating differential expression, this time keeping all original
#proteins (for volcano plot)
#and changing p.value cut-off to 0.10
calcDiffExpr2<-function(label.g1, label.g2){
  message("calcDiffExpr: ",label.g1,"|",label.g2)
  g1<-prot[,grepl(label.g1,colnames(prot))]
  g2<-prot[,grepl(label.g2,colnames(prot))]
  
  prot.sel<-cbind(g1, g2)
  group<-c(rep(label.g1, ncol(g1)), rep(label.g2, ncol(g2)))
  
  res_anova<-apply(prot.sel, 1, function(x){
    pval=doANOVA(group=group, data.row=x)
    return(data.frame(pval=pval))
  })
  
  res_logfc <- apply(prot.sel, 1, function(x){
    mean.g1=mean(x[group==unique(group)[1]])
    mean.g2=mean(x[group==unique(group)[2]])
    fold.change<-(2^mean.g1-2^mean.g2) / 2^mean.g1 # expression values in log scale
    return(data.frame(fold.change=fold.change))
  })
  
  df_anova <- data.frame(gene=rownames(prot.sel), bind_rows(res_anova))
  df_logfc <- data.frame(gene=rownames(prot.sel), bind_rows(res_logfc))
  # 
  df <- cbind(df_anova, df_logfc)
  
  #get rid of duplicate gene column
  df[,3]<-NULL
  
  #Bonferroni multiple testing adjustment
  #multiplies the raw P values by the number of tests m (i.e. length of the vector P_values). Using the p.adjust 
  #function and the 'method' argument set to "bonferroni", we get a vector of same length but with adjusted P values.
  #This adjustment approach corrects according to the family-wise error rate of at least one false positive
  #(FamilywiseErrorRate (FWER)=Probability (FalsePositive???1))
  
  df$pval.adj.bonferroni <- p.adjust(df$pval, method="bonferroni")
  
  df.sel <- df[df$pval<0.1 & abs(df$fold.change)>0.5,]
  
  #get rid of rows where the gene column is NA
  df.sel_no_NAs <- df.sel %>% drop_na(gene)
  
  ####adjust p.values to correct for multiple testing
  
  #Benjamini-Hochberg/FDR
  #false discovery rate (FalseDiscoveryRate (FDR)=Expected (FalsePositive/ (FalsePositive+TruePositive))). 
  #In other words, FDR is the expected proportion of false positives among all positives which rejected 
  #the null hypothesis and not among all the tests undertaken. In the FDR method, P values are ranked in an 
  #ascending array and multiplied by m/k where k is the position of a P value in the sorted vector and m is 
  #the number of independent tests. 
  
  df.sel_no_NAs$pval.adj.fdr <- p.adjust(df.sel_no_NAs$pval, method="fdr")
  df.sel_no_NAs <- df.sel_no_NAs
  
  #fname.out<-file.path(dir.results, "04_diff_gene_expression",paste0("diff_",label.g1,"_",label.g2,".csv"))
  
  #remove the genes which have "-Inf" in the fold.change column
  df_final <-df.sel_no_NAs[!(df.sel_no_NAs$fold.change == "-Inf"),]
  
  saveRDS(df_final, file = paste0("diff_expr_padj_0.10_all",label.g1,"_",label.g2,".rds"))
  fname.out<-file.path(dir.results, paste0("diff_expr_padj_0.10_all",label.g1,"_",label.g2,".csv"))
  
  write.csv(df_final, file=fname.out,row.names=FALSE, quote = FALSE)
}


#function for plotting ANOVA and log 2 Fold Change values on a volcano plot

volcano_ggplot <- function(top_table, p_value_cutoff, subtitle){
  ggplot(data=top_table, aes(x=fold.change, y=-log10(pval.adj.fdr), label=gene)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +    
    scale_color_manual(values=c("#F68D91",  "#4E80A5", "grey", "#A4D49C", "orange")) +
    geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(p_value_cutoff), col="red") +
    labs(title = 'Organoid Protein Secretome dataset',
         subtitle=paste0(subtitle, " p. adj. value cut-off = ", p_value_cutoff),
         x = 'log2 Fold Change')
}
