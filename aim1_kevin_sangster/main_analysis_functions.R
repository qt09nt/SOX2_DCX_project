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

volcano_plot = function(top_table, title){
  v = EnhancedVolcano(top_table, x="logFC", y= "adj.P.Val", lab=row.names(top_table), title = title,
                      pCutoff = 0.1, FCcutoff = 1.333, 
                      ylim = c(0, max(-log10(top_table$adj.P.Val), na.rm=TRUE)),
                      col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
                      colAlpha = 1.0)
  pdf(paste("figures/", title, "volcano.pdf"), width = 8, height = 8)
  print(v)
  dev.off()
  return(v)
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




