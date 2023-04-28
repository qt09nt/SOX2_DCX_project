#functions for analysis of single cell RNA seq data from Ramos et al (2022)

#plot the differentially expressed genes on a volcano plot
volcano_ggplot <- function(top_table, p_value_cutoff, subtitle){
  ggplot(data=top_table, aes(x=avg_log2FC, y=-log10(p_val_adj), label=gene)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +    
    scale_color_manual(values=c("#F68D91",  "#4E80A5", "grey", "#A4D49C", "orange")) +
    geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(p_value_cutoff), col="red") +
    labs(title = 'Ramos et al (2022) GSM6720853 Glutamatergic Neurons with high FABP7 vs Glutamatergic Neurons with low FABP7',
         subtitle=paste0(subtitle, " p. adj. value cut-off = ", p_value_cutoff),
         x = 'log2 Fold Change')
}