#functions used to analyze BrainSpan_36_donors_TPM.R

#facet boxplots of 
#re-write RNA Seq facet boxplots of ganglionic eminence vs cortex at 8pcw and 9pcw as a function so there there's less code repetition:
plotting <- function(dataframe, gene){
  gene_name <- toString(gene)
  gene <- dataframe %>%    #extract the relevant column of RNA seq TPM values for that gene
    pull(gene)
  r <<- ggplot(data = dataframe, aes(x = region,  y = gene, group = region, colour = region)) +
    geom_boxplot() +
    geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
    facet_wrap(~ age) +  #separate the boxplots by age
    theme_bw() + #set the background to white
    theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    labs(title = 'BrainSpan Transcriptome Dataset',
         subtitle=paste0(gene_name, ' expression in ganglionic eminence and cortex at 8pcw and 9pcw'),
         x = 'structure',
         y = paste0(gene_name,' expression (TPM)'))
  
  #Adding p value: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  #  Add p-value
  # Change method
  r + stat_compare_means(method = "t.test", aes(hjust = 50, vjust = 70)) + scale_color_manual(values=c("#F69679", "#9DB4DD"))
}


plotting_cortex <- function(dataframe, gene){
  gene_name <- toString(gene)
  gene <- dataframe %>%    #extract the relevant column of RNA seq TPM values for that gene
    pull(gene)
  ggplot(data = cortex, aes(x = structure_name,  y = gene, colour = structure_name)) +
    geom_boxplot() +
    geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
    facet_wrap(~ age) +  #separate the boxplots by age
    theme_bw() + #set the background to white
    theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title = 'BrainSpan Transcriptome Dataset',
         subtitle=paste0(gene_name, ' expression in cortex at 8pcw to 40 years'),
         x = 'structure',
         y = paste0(gene_name,' expression (TPM)'))
}


#dataframe should be combined RNA expression in TPM values combined with columns meta and filtered for 
#ganglionic eminences and cortex
plotting_ge_cortex <- function(dataframe, gene){
  gene_name <- toString(gene)
  gene <- dataframe %>%    #extract the relevant column of RNA seq TPM values for that gene
    pull(gene)
  ggplot(data = dataframe, aes(x = structure_name,  y = gene, colour = structure_name)) +
    geom_boxplot() +
    geom_jitter(width = 0.10) +      #show the separate data points using jitter to slightly separate the points to visualize easier
    facet_wrap(~ age) +  #separate the boxplots by age
    theme_bw() + #set the background to white
    theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title = 'BrainSpan Transcriptome Dataset',
         subtitle=paste0(gene_name, ' expression in all brain structures at 8pcw to 40 years'),
         #subtitle=paste0(gene_name, ' expression in cortex and ganglionic eminences at 8pcw and 9pcw'),
         x = 'structure',
         y = paste0(gene_name,' expression (TPM)'))
}

#function for extracting particular time points for comparison of FABP7 expression in prenatal brain (BrainSpan data)
#input "combined_df" is the dataframe "combined" which contains the FABP7 expression values in one column
#as well as the sample meta data ie. brain structure, age etc
get_tp <-function(combined_df, timepoint1, timepoint2){
  tp1 <- combined_df[combined_df$age == timepoint1, ]
  tp2 <- combined_df[combined_df$age == timepoint2, ]
  both <- rbind(tp1, tp2)
  return(both)
}

#function for imputing FABP7 expression values which are missing from certain structures in BRainSPan
#test_df is the dataframe forFABP7 expression containing 2 timepoints to be tested
impute_process <-function(test_df){
  #impute NA values
  ini <- mice(data = test_df, maxit = 0)
  
  #define the 'subject' identifier (code as '-2' in predictor matrix)
  pred <- ini$pred
  pred["FABP7", "structure_ID"] <- -2
  
  #run MI
  #"2l.pan is level-1 normal homoscedastic, lmer
  imp <- mice(data = test_df, pred = pred, method = "2l.pan", m = 100)
  
  summary(imp)
  
  #create a list of completed data sets
  implist <<- mids2mitml.list(imp)
  
  #fit the ANOVA model using lmer(). The model contains an additional term (1|electrode_number), which
  #specifies a random effect for each subject. This effect captures unsystematic differences between subjects,
  #thus accounting for the nested structure of the repeated-measures data.
  
  #fit the ANOVA model
  fit3 <<- with(implist, lmer(FABP7 ~ 1 + age + (1|structure_ID)))
  testEstimates(fit3, var.comp = TRUE)
}

#fit ANOVA model with lmer (fit reduced model)
fit_anova<-function(implist){
  #pool the parameter estimates
  fit3.reduced <- with(implist, lmer(FABP7 ~ 1 + (1|structure_ID)))
  testModels(fit3, fit3.reduced, method = "D1")
}
