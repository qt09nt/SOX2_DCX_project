#July 6, 2022

#DCX timepoint 1 comparison to DCX timepoint 3
#volcano plot and pathway analysis

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/")

source("R/perseus_functions.R")
source("R/main_analysis_functions.R")
source("R/protein_organoid_working_matrix_main_analysis_functions.R")

library(tidyverse)
library(viridis)
library(forcats)    #used for handling factors for gsea pathways dataframe & reordering by FDR.qvalue
library(dplyr)
library(stringr)  #for subsetting gsea results for separate Gene Ontology Molecular Functions (GOMF), Gene Ontology Cellular Components (GOCC)
                    #and GOBP (Biological Processes) plots
library(pheatmap)
library(MCMCglmm)
library(ggplot2)
library(EnhancedVolcano)
library(limma)
library(factoextra)
library(dendextend)
library(purrr)
library(gridExtra)
library(cowplot)
library(hrbrthemes)
library(viridis)
library(WGCNA)
library(RColorBrewer)
library(ggsignif)
library(colorspace)
library(ggpubr)
library(cluster)
library(factoextra)
library(DEGreport)

#################

#read in previously processed data (ie. where only YFP+ samples are kept)
data_pos <- read.table(file = "results/data_pos.txt", sep = "\t", header = TRUE, row.names = 1)
metadata_pos <- read.table(file = "results/metadata_pos.txt", sep = "\t", header = TRUE)

#rename the colnames of data_pos (YFP +) to match the sample names in metadata column
colnames(data_pos) <- (metadata_pos$sample)

#first filter for DCX samples only in the metadata table, and DCX YFP+ in protein expression matrix:
DCX_meta <- metadata_pos[metadata_pos$cell_type == "DCX",]
DCX_expr <- data_pos[,1:7]

#keep only DCX time point 1 & DCX time point 3 in metadata table
DCX_1_3_meta <- DCX_meta[(DCX_meta$tp == 1 |DCX_meta$tp == 3),]

#filter protein data matrix for DCX timepoints 1 & 3 only:
dcx_tp_1_3_samples <-c(DCX_1_3_meta$sample)

dcx_tp_1_3_expr <- DCX_expr %>% 
  select(all_of(dcx_tp_1_3_samples))

#remove redundant unfiltered metadata and expression matrices containing
#all DCX samples from env to prevent confusion
rm(DCX_meta)
rm(DCX_expr)
rm(metadata_pos)

#do differential expression analysis between DCX YFP+ timepoint 1 and DCX YFP+ timepoint 3, using limma package
difference = ifelse(DCX_1_3_meta$tp == "1", 1, 0)
design = cbind(control = 1, difference = difference)
dcx_tp_1_vs_3_top_table = get_top_table(dcx_tp_1_3_expr, design)

#save diff expression top table results
write.csv(dcx_tp_1_vs_3_top_table, "results/dcx_tp_1_vs_3_top_table.csv")

#adjust FDR cutoffs in 
#C:\Users\qt09n\Desktop\Technical Analyst I UHN May 4 2021\organoid group\Sofia\pathway analysis\R\main_analysis_functions.R
#adjusted pvalue cutoff 0.1:
print(volcano_plot(dcx_tp_1_vs_3_top_table, "July 8 2022 Organoid DCX Samples YFP+ Time Point 1 vs 3 adj pvalue cutoff 0.1"))

#adj.pvalue cutoff 0.01
print(volcano_plot(dcx_tp_1_vs_3_top_table, "July 8 2022 Organoid DCX Samples YFP+ Time Point 1 vs 3, adj.pvalue cutoff 0.01"))

#adj. pvalue cutoff 0.05
print(volcano_plot(dcx_tp_1_vs_3_top_table, "July 8 2022 Organoid DCX YFP+ TimePoint 1 vs TimePoint 3, adj.pvalue cutoff 0.05"))

#adj. pvalue cutoff 0.25
print(volcano_plot(dcx_tp_1_vs_3_top_table, "July 8 2022 Organoid DCX YFP+ TimePoint 1 vs TimePoint 3, adj.pvalue cutoff 0.25"))



###### re-do volcano plots using ggplot to better visualize the significant proteins

#first label which proteins are differentially expressed according to adj.pvalue cutoffs & log2FC cutoffs:
dcx_tp_1_vs_3_top_table_labeled <- get_diffexpressed_modified(top_table = dcx_tp_1_vs_3_top_table, 
                                                              p_value=0.1, condition1 ="DCX timepoint 1", condition2 = "DCX timepoint 3")


#label genes which are significantly up in DCX tp 1 and significantly up in DCX tp 3:
dcx_tp_1_vs_3_top_table_labeled$genelabels <- ifelse(dcx_tp_1_vs_3_top_table_labeled$diffexpressed == "UP in DCX timepoint 1" 
                                                    | dcx_tp_1_vs_3_top_table_labeled$diffexpressed == "UP in DCX timepoint 3"
                                                   , TRUE, FALSE)
  
modified_volcano_ggplot(top_table=dcx_tp_1_vs_3_top_table_labeled, 0.10, "DCX timepoint 1 vs DCX timepoint 3")


write.csv(dcx_tp_1_3_expr, "results/dcx_tp_1_3_expr.csv")









#########################double check with Sofia's Perseus differential expression analysis results
perseus_diffexpr <- read.table(file = "results/perseus_DCX_TP1_VS_TP3.txt", sep = "\t", header = TRUE, fill=TRUE)

#filter perseus diff expr results for significant genes only:
perseus_significant <- perseus_diffexpr[perseus_diffexpr$Significant == "+",]

#export table for significant proteins from Perseus results:
write.csv(perseus_significant, "results/perseus_significant_dcx_tp1_vs_tp3.csv")


#check if significant proteins found using Limma are significant in volcano plots from Perseus
any(perseus_significant$Gene.names=="COX7A2")  #False: NOT FOUND In perseus as significant
any(perseus_significant$Gene.names=="LGALS3")  #False: NOT FOUND In perseus as significant
any(perseus_significant$Gene.names=="LMAN1")   #False: NOT FOUND In perseus as significant
any(perseus_significant$Gene.names=="FDFT1")  #False: NOT FOUND In perseus as significant
any(perseus_significant$Gene.names=="CCAR1")  #False: NOT FOUND In perseus as significant




#################### double check with Perseus protein matrix that Sofia used to generate volcano plot




############ DCX Time Point 1 vs DCX Time Point 2:
#use get_timepoints function

#adj.p.value cutoff = 0.1
get_timepoints(annotations = DCX_meta, working_data = DCX_expr, timepoint1 = 1, timepoint2 = 2)

#adj.p.value cutoff = 0.25
get_timepoints(annotations = DCX_meta, working_data = DCX_expr, timepoint1 = 1, timepoint2 = 2)

###### re-do volcano plots using ggplot to better visualize the significant proteins

#read in DCX tp 1 vx tp 2 top table
dcx_tp_1_vs_2_top_table <- read.table(file = "results/DCX1_vs_2_top_table.txt", sep = "\t", header = TRUE, row.names = 1)


#first label which proteins are differentially expressed according to adj.pvalue cutoffs & log2FC cutoffs:
dcx_tp_1_vs_2_top_table_labeled <- get_diffexpressed_modified(top_table = dcx_tp_1_vs_2_top_table, 
                                                              p_value=0.25, condition1 ="DCX timepoint 1", condition2 = "DCX timepoint 2")


#label genes which are significantly up in DCX tp 1 and significantly up in DCX tp 2:
dcx_tp_1_vs_2_top_table_labeled$genelabels <- ifelse(dcx_tp_1_vs_2_top_table_labeled$diffexpressed == "UP in DCX timepoint 1" 
                                                     | dcx_tp_1_vs_2_top_table_labeled$diffexpressed == "UP in DCX timepoint 2"
                                                     , TRUE, FALSE)

dev.off()
dev.new()
modified_volcano_ggplot(top_table=dcx_tp_1_vs_2_top_table_labeled, 0.25, "DCX timepoint 1 vs DCX timepoint 2")


#load in GSEA results for DCX timepoint 1 vs DCX timepoint 2
#metric used for ranking genes = "Diff of Classes"
filter_gsea_DCX_tps(gsea_report_na_pos = "results/dcx_tp2_vs_tp1_GO_human_all_pathways.Gsea.1657654730497/gsea_report_for_tp1_1657654730497.tsv",
                    gsea_report_for_na_neg ="results/dcx_tp2_vs_tp1_GO_human_all_pathways.Gsea.1657654730497/gsea_report_for_tp2_1657654730497.tsv", timepoint1="timepoint1", timepoint2="timepoint2")
  
#filter for GOBP database results
DCX_tp1_GOBP <- filter(filtered_DCX, grepl("GOBP", filtered_DCX$Pathways))
DCX_tp2_GOBP <- filter(filtered_DCX2, grepl("GOBP", filtered_DCX2$Pathways))

#94 and 201 pathways respectively for DCX tp 1 & DCX tp 2 at FDR.q.value of 0.10

#use more stringent FDR to further filter pathways -> try FDR.q.value of 0.05
DCX_tp1_GOBP2<- DCX_tp1_GOBP[DCX_tp1_GOBP$FDR.q.val <= 0.05, ]

#get rid of 3rd column which is empty:
DCX_tp1_GOBP2[,3]<- NULL


####### subset protein matrix for DCX timepoint 1 & DCX timepoint 2 for GSEA analysis:

#keep only DCX time point 2 & DCX time point 3 in metadata table
DCX_1_2_meta <- DCX_meta[(DCX_meta$tp == 1 |DCX_meta$tp == 2),]

#filter protein data matrix for DCX timepoints 1 & 2 only:
dcx_tp_1_2_samples <-c(DCX_1_2_meta$sample)

dcx_tp_1_2_expr <- DCX_expr %>% 
  select(all_of(dcx_tp_1_2_samples))

write.csv(dcx_tp_1_2_expr, "results/dcx_tp_1_2_expr.csv", row.names = TRUE)




################################## DCX TimePoint 2 vs DCX Time Point 3:

#adj.p.value cutoff = 0.25
get_timepoints(annotations = DCX_meta, working_data = DCX_expr, timepoint1 = 2, timepoint2 = 3)

#adj.p.value cutoff = 0.1
get_timepoints(annotations = DCX_meta, working_data = DCX_expr, timepoint1 = 2, timepoint2 = 3)

###### re-do volcano plots using ggplot to better visualize the significant proteins

#read in DCX tp 2 vx tp 3 top table
dcx_tp_2_vs_3_top_table <- read.table(file = "results/DCX tp2_vs_tp3_top_table.txt", sep = "\t", header = TRUE, row.names = 1)


#first label which proteins are differentially expressed according to adj.pvalue cutoffs & log2FC cutoffs:
dcx_tp_2_vs_3_top_table_labeled <- get_diffexpressed_modified(top_table = dcx_tp_2_vs_3_top_table, 
                                                              p_value=0.25, condition1 ="DCX timepoint 2", condition2 = "DCX timepoint 3")


#label genes which are significantly up in DCX tp 1 and significantly up in DCX tp 3:
dcx_tp_2_vs_3_top_table_labeled$genelabels <- ifelse(dcx_tp_2_vs_3_top_table_labeled$diffexpressed == "UP in DCX timepoint 2" 
                                                     | dcx_tp_2_vs_3_top_table_labeled$diffexpressed == "UP in DCX timepoint 3"
                                                     , TRUE, FALSE)

modified_volcano_ggplot(top_table=dcx_tp_2_vs_3_top_table_labeled, 0.25, "DCX timepoint 2 vs DCX timepoint 3")


#########get expression data for DCX tp 2 and DCX tp 3 for GSEA analysis:

#keep only DCX time point 2 & DCX time point 3 in metadata table
DCX_2_3_meta <- DCX_meta[(DCX_meta$tp == 2 |DCX_meta$tp == 3),]

#filter protein data matrix for DCX timepoints 2 & 3 only:
dcx_tp_2_3_samples <-c(DCX_2_3_meta$sample)

dcx_tp_2_3_expr <- DCX_expr %>% 
  select(all_of(dcx_tp_2_3_samples))

write.csv(dcx_tp_2_3_expr, "results/dcx_tp_2_3_expr.csv", row.names = TRUE)


############### pathway analysis using GSEA, databases used: GOCC, GOMF, GOBP and Human All Pathways
#not enough samples in DCX time point 2 (2 replicates only), and time point 3 (2 replicates only)
#to do tTest as metric for ranking, so did Diff_of_classes as ranking metric
#load in gsea results:

###########DCX tp 3 vs DCX tp 1

#DCX tp 1 set nom.pvalue and fdr.q.value to 0.20
#DCX tp 3 set nom.pvalue and fdr.q.value to 0.05
filter_gsea_DCX_tps(gsea_report_na_pos="results/dcx_tp3_vs_tp1_GO_human_all_pathways.Gsea.1657653244475/gsea_report_for_tp1_1657653244475.tsv", 
                                gsea_report_for_na_neg="results/dcx_tp3_vs_tp1_GO_human_all_pathways.Gsea.1657653244475/gsea_report_for_tp3_1657653244475.tsv",
                    timepoint1 = "timepoint 1", timepoint2="timepoint 3")

#filter for GOCC results only  
extract_specific_GO_DCX (filtered_DCX=filtered_DCX, filtered_DCX2=filtered_DCX2,
                         GO_db="GOCC", timepoint1="timepoint 1", timepoint2="timepoint 3")

#plot results
ggplot(GO_merged, aes(x = timepoint, y = Pathways, size = NES, color = FDR.q.val )) +
  geom_point(alpha = 0.7) + #scatterplot
  scale_size(range = c(2, 5), name = "NES") +  #change the size of the points and create a name for the size legend
  scale_color_viridis(discrete = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle("GO Cellular Components GSEA Pathway Analysis - DCX timepoint 1 vs DCX timepoint 3")

#filter for GOBP:
extract_specific_GO_DCX (filtered_DCX=filtered_DCX, filtered_DCX2=filtered_DCX2,
                         GO_db="GOBP", timepoint1="timepoint 1", timepoint2="timepoint 3")

#remove some of the redundant pathways to make it easier to read on the plot y-axis
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_ACTIN_FILAMENT_BASED_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_EPITHELIAL_CELL_DEVELOPMENT",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_LOCOMOTION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_CELL_JUNCTION_ASSEMBLY", ]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_WOUND_HEALING",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_MRNA_PROCESSING",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_RNA_SPLICING",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_ACTIN_FILAMENT_BUNDLE_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_RNA_PROCESSING",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_MRNA_METABOLIC_PROCESS",]

#group pathways together by timepoint
tp1 <- GO_merged[GO_merged$timepoint == "timepoint 1", ]
tp1$Pathways

tp3 <- GO_merged[GO_merged$timepoint == "timepoint 3", ]
tp3$Pathways

GO_merged$Pathways <- factor(GO_merged$Pathways, levels = c( "GOBP_MRNA_PROCESSING",                                                        
 "GOBP_RNA_SPLICING_VIA_TRANSESTERIFICATION_REACTIONS",                         
 "GOBP_MRNA_METABOLIC_PROCESS",                                                 
 "GOBP_MRNA_3_END_PROCESSING",                                                  
 "GOBP_RNA_POLYADENYLATION",                                                    
 "GOBP_NEGATIVE_REGULATION_OF_MRNA_METABOLIC_PROCESS",                          
 "GOBP_SPLICEOSOMAL_COMPLEX_ASSEMBLY",                                          
 "GOBP_POSITIVE_REGULATION_OF_RNA_SPLICING",                                    
 "GOBP_PROTEIN_DNA_COMPLEX_SUBUNIT_ORGANIZATION",                               
 "GOBP_NEGATIVE_REGULATION_OF_MRNA_CATABOLIC_PROCESS",                          
 "GOBP_TRANSCRIPTION_ELONGATION_FROM_RNA_POLYMERASE_II_PROMOTER",               
 "GOBP_RNA_3_END_PROCESSING",                                                   
 "GOBP_NUCLEOSOME_ORGANIZATION",                                                
 "GOBP_CHROMATIN_ORGANIZATION",                                                 
 "GOBP_PRODUCTION_OF_MOLECULAR_MEDIATOR_OF_IMMUNE_RESPONSE",                    
 "GOBP_POSITIVE_REGULATION_OF_MRNA_PROCESSING",                                 
 "GOBP_MRNA_CIS_SPLICING_VIA_SPLICEOSOME",                                      
 "GOBP_NEGATIVE_REGULATION_OF_NUCLEOBASE_CONTAINING_COMPOUND_METABOLIC_PROCESS",
 "GOBP_ISOPRENOID_METABOLIC_PROCESS",                                           
 "GOBP_NUCLEAR_TRANSCRIBED_MRNA_CATABOLIC_PROCESS",                             
 "GOBP_DEFENSE_RESPONSE_TO_SYMBIONT",                                           
 "GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_GROWTH_FACTOR_STIMULUS",
 "GOBP_BLOOD_VESSEL_MORPHOGENESIS",                              
  "GOBP_TUBE_MORPHOGENESIS",                                    
 "GOBP_BIOLOGICAL_ADHESION",                                     
 "GOBP_VASCULATURE_DEVELOPMENT",                                 
 "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION",           
  "GOBP_TUBE_DEVELOPMENT",                                        
  "GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION",                    
  "GOBP_CELL_SUBSTRATE_ADHESION",                                 
 "GOBP_ANATOMICAL_STRUCTURE_FORMATION_INVOLVED_IN_MORPHOGENESIS",
 "GOBP_RESPONSE_TO_WOUNDING",                                    
 "GOBP_CELL_CELL_ADHESION",                                      
 "GOBP_CIRCULATORY_SYSTEM_PROCESS",                              
  "GOBP_MONONUCLEAR_CELL_MIGRATION",                              
  "GOBP_GLIOGENESIS",                                             
  "GOBP_REGULATION_OF_CELL_ADHESION",                             
 "GOBP_GLIAL_CELL_DIFFERENTIATION",                              
 "GOBP_VASCULAR_PROCESS_IN_CIRCULATORY_SYSTEM",                  
 "GOBP_ACTOMYOSIN_STRUCTURE_ORGANIZATION",                       
 "GOBP_CELL_MIGRATION",                                         
 "GOBP_REGULATION_OF_CELLULAR_COMPONENT_MOVEMENT",               
 "GOBP_CELL_CHEMOTAXIS",                                         
 "GOBP_ENDOTHELIAL_CELL_MIGRATION",                              
  "GOBP_CONNECTIVE_TISSUE_DEVELOPMENT",                           
  "GOBP_ACTIN_FILAMENT_ORGANIZATION",                             
 "GOBP_POSITIVE_REGULATION_OF_CELL_SUBSTRATE_ADHESION",          
 "GOBP_EPITHELIAL_CELL_DIFFERENTIATION",                         
 "GOBP_REGULATION_OF_LEUKOCYTE_MIGRATION",                       
 "GOBP_SUPRAMOLECULAR_FIBER_ORGANIZATION",                       
 "GOBP_CIRCULATORY_SYSTEM_DEVELOPMENT",                          
 "GOBP_STRIATED_MUSCLE_CELL_DIFFERENTIATION",
 "GOBP_EPITHELIUM_DEVELOPMENT",                                  
 "GOBP_TISSUE_DEVELOPMENT",                                      
 "GOBP_CELL_MATRIX_ADHESION",                                    
 "GOBP_LEUKOCYTE_MIGRATION",                                     
 "GOBP_POSITIVE_REGULATION_OF_LOCOMOTION",                       
 "GOBP_METAL_ION_TRANSPORT",                                     
 "GOBP_AMEBOIDAL_TYPE_CELL_MIGRATION",                           
 "GOBP_REGULATION_OF_ENDOTHELIAL_CELL_MIGRATION",                
 "GOBP_CELL_CELL_JUNCTION_ORGANIZATION",                         
 "GOBP_PROTEIN_LOCALIZATION_TO_PLASMA_MEMBRANE",                 
 "GOBP_PHAGOCYTOSIS",                                            
 "GOBP_PROTEIN_LOCALIZATION_TO_CELL_PERIPHERY",                  
 "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",                     
 "GOBP_ANION_TRANSMEMBRANE_TRANSPORT",                           
 "GOBP_SENSORY_PERCEPTION_OF_MECHANICAL_STIMULUS",               
 "GOBP_MONONUCLEAR_CELL_DIFFERENTIATION"))

plot_NES_merged(filtered_df=GO_merged, GO_db="GOBP")

############ DCX timepoint 1 vs DCX timepoint 3 GSEA GOMF


extract_specific_GO_DCX (filtered_DCX=filtered_DCX, filtered_DCX2=filtered_DCX2,
                         GO_db="GOMF", timepoint1="timepoint 1", timepoint2="timepoint 3")

plot_NES_merged(GO_merged, "GOMF")

#try with less stringent FDR cutoff for GOMF 
#FDR.q.val` <= 0.25 for DCX timepoint 1
#FDR.q.val` <= 0.1 for DCX timepoint 3
filter_gsea_DCX_tps(gsea_report_na_pos="results/dcx_tp3_vs_tp1_GO_human_all_pathways.Gsea.1657653244475/gsea_report_for_tp1_1657653244475.tsv", 
                    gsea_report_for_na_neg="results/dcx_tp3_vs_tp1_GO_human_all_pathways.Gsea.1657653244475/gsea_report_for_tp3_1657653244475.tsv",
                    timepoint1 = "timepoint 1", timepoint2="timepoint 3")



############################## Pathway analysis for DCX TimePoint 2 vs DCX Timepoint 3

#load in GSEA results

#DCX tp 2 set nom.pvalue and fdr.q.value to  & NOM p.value to 0.25
#DCX tp 3 set nom.pvalue and fdr.q.value to 0.25 & NOM pvalue to 0.35
filter_gsea_DCX_tps(gsea_report_na_pos="results/dcx_tp2_vs_tp3_GO_human_all_pathways.Gsea.1657652459317/gsea_report_for_tp2_1657652459317.tsv", 
                    gsea_report_for_na_neg="results/dcx_tp2_vs_tp3_GO_human_all_pathways.Gsea.1657652459317/gsea_report_for_tp3_1657652459317.tsv",
                    timepoint1 = "timepoint 2", timepoint2="timepoint 3")

#filter for GOCC:
extract_specific_GO_DCX (filtered_DCX=filtered_DCX, filtered_DCX2=filtered_DCX2,
                         GO_db="GOCC", timepoint1="timepoint 2", timepoint2="timepoint 3")

#group all the timepoint 2 pathways together to make it easier to read from plot
DCX_GO$Pathways
DCX2_GO$Pathways

GO_merged$Pathways <- factor(GO_merged$Pathways, levels = c("GOCC_AXON",                                                      
 "GOCC_MICROTUBULE_ASSOCIATED_COMPLEX",                            
 "GOCC_MICROTUBULE",                                               
 "GOCC_DISTAL_AXON",                                               
 "GOCC_NEURON_PROJECTION",                                         
 "GOCC_INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE",                    
 "GOCC_SITE_OF_POLARIZED_GROWTH",                                  
 "GOCC_SOMATODENDRITIC_COMPARTMENT",                               
 "GOCC_NEURON_TO_NEURON_SYNAPSE",                                  
 "GOCC_CELL_BODY",                                                 
 "GOCC_PLASMA_MEMBRANE_REGION",                                    
 "GOCC_CYTOPLASMIC_REGION",                                        
 "GOCC_POLYMERIC_CYTOSKELETAL_FIBER",                              
 "GOCC_RECYCLING_ENDOSOME",                                        
 "GOCC_PRESYNAPSE",                                                
 "GOCC_SYNAPSE",                                                   
 "GOCC_CILIUM",                                                    
 "GOCC_DENDRITIC_TREE",                                            
 "GOCC_ACTIN_BASED_CELL_PROJECTION",                               
 "GOCC_APICAL_PLASMA_MEMBRANE",                                    
 "GOCC_INCLUSION_BODY",                                            
 "GOCC_TRANS_GOLGI_NETWORK",                                       
 "GOCC_POSTSYNAPSE",                                               
 "GOCC_NEURON_SPINE",                                              
 "GOCC_EXTRINSIC_COMPONENT_OF_CYTOPLASMIC_SIDE_OF_PLASMA_MEMBRANE",
 "GOCC_FILOPODIUM",                                                
 "GOCC_BASOLATERAL_PLASMA_MEMBRANE",                               
 "GOCC_SUPRAMOLECULAR_POLYMER",                                    
 "GOCC_APICAL_PART_OF_CELL",                                       
 "GOCC_TRANSPORTER_COMPLEX",                                       
 "GOCC_NEURON_PROJECTION_CYTOPLASM",                               
 "GOCC_CELL_PROJECTION_MEMBRANE",                                  
 "GOCC_LEADING_EDGE_MEMBRANE",                                     
 "GOCC_AXON_CYTOPLASM",                                            
 "GOCC_EXOCYTIC_VESICLE",                                          
  "GOCC_PLASMA_MEMBRANE_PROTEIN_COMPLEX",
 "GOCC_NUCLEAR_CHROMOSOME",
 "GOCC_U2_TYPE_SPLICEOSOMAL_COMPLEX",      
 "GOCC_NUCLEAR_PROTEIN_CONTAINING_COMPLEX", "GOCC_NUCLEOLUS",                         
 "GOCC_LYSOSOMAL_LUMEN",                    "GOCC_PRECATALYTIC_SPLICEOSOME",          
 "GOCC_SPLICEOSOMAL_COMPLEX",               "GOCC_FIBRILLAR_CENTER",                  
 "GOCC_NUCLEAR_PERIPHERY",                  "GOCC_CYTOPLASMIC_STRESS_GRANULE",       
 "GOCC_METHYLTRANSFERASE_COMPLEX",          "GOCC_SPLICEOSOMAL_SNRNP_COMPLEX",        
 "GOCC_P_BODY",                             "GOCC_NUCLEAR_SPECK",                     
 "GOCC_POSTSYNAPTIC_MEMBRANE",              "GOCC_RIBONUCLEOPROTEIN_COMPLEX",         
 "GOCC_SM_LIKE_PROTEIN_FAMILY_COMPLEX"))   

#plot merged GOCC pathways for DCX timepoint 2 & DCX timepoint 3
plot_NES_merged(GO_merged, GO_db ="GOCC", firsttimepoint = "2", secondtimepoint = "3")

###########filter for GOMF for DCX tp 2 vs DCX tp 3 comparison:
extract_specific_GO_DCX (filtered_DCX=filtered_DCX, filtered_DCX2=filtered_DCX2,
                         GO_db="GOMF", timepoint1="timepoint 2", timepoint2="timepoint 3")

#group tp2 pathways together and tp3 pathways together
DCX_GO$Pathways
DCX2_GO$Pathways

GO_merged$Pathways <- factor(GO_merged$Pathways, levels = c(
  "GOMF_GDP_BINDING",                                                   "GOMF_GTPASE_ACTIVITY",                                              
  "GOMF_ACTIVE_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",                 "GOMF_CYTOSKELETAL_MOTOR_ACTIVITY",                                  
  "GOMF_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",                     "GOMF_MICROTUBULE_BINDING",                                          
  "GOMF_CATION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",                     "GOMF_INORGANIC_MOLECULAR_ENTITY_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
  "GOMF_PRIMARY_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",             "GOMF_ATPASE_BINDING",                                               
  "GOMF_GUANYL_NUCLEOTIDE_BINDING",                                     "GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",                       
  "GOMF_HEXOSYLTRANSFERASE_ACTIVITY",                                   "GOMF_TUBULIN_BINDING", 
  "GOMF_DNA_HELICASE_ACTIVITY",                       "GOMF_ATP_DEPENDENT_ACTIVITY_ACTING_ON_DNA",       
  "GOMF_HELICASE_ACTIVITY",                           "GOMF_CATALYTIC_ACTIVITY_ACTING_ON_A_NUCLEIC_ACID",
  "GOMF_CATALYTIC_ACTIVITY_ACTING_ON_DNA",            "GOMF_NUCLEASE_ACTIVITY",                          
  "GOMF_NUCLEOSOMAL_DNA_BINDING",                     "GOMF_NUCLEOSOME_BINDING")) 

plot_NES_merged(GO_merged, GO_db ="GOMF", firsttimepoint = "2", secondtimepoint = "3")


######### DCX tp 2 vs DCX tp 3 pathway analysis in GOBP

extract_specific_GO_DCX (filtered_DCX=filtered_DCX, filtered_DCX2=filtered_DCX2,
                         GO_db="GOBP", timepoint1="timepoint 2", timepoint2="timepoint 3")

#remove some redundant pathways
#can use a word frequency counter to help see which key words/pathways are frequent and remove redundant ones:
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_CELL_JUNCTION_ASSEMBLY",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_AXON_EXTENSION",]

GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CELL_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CELL_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_CELL_PROJECTION_ASSEMBLY",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_MICROTUBULE_BASED_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_MICROTUBULE_POLYMERIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_MICROTUBULE_CYTOSKELETON_ORGANIZATION", ]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_HEART_CONTRACTION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_HEART_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_REGULATION_OF_CELL_PROJECTION_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_EXTENT_OF_CELL_GROWTH",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_REGULATION_OF_DEVELOPMENTAL_GROWTH",]



GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_REGULATION_OF_ION_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_CATION_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_CELL_PROJECTION_ASSEMBLY",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CELL_PROJECTION_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_PROTEIN_POLYMERIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_REGULATION_OF_MICROTUBULE_POLYMERIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_CELL_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_MICROTUBULE_CYTOSKELETON_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CELL_PROJECTION_ASSEMBLY",]

GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_INORGANIC_ION_TRANSMEMBRANE_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_NEURON_DIFFERENTIATION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_DEVELOPMENTAL_CELL_GROWTH",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_POST_GOLGI_VESICLE_MEDIATED_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CELL_PROJECTION_ASSEMBLY",]

GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_CELL_PROJECTION_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_CELL_MORPHOGENESIS",]


GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_CELL_JUNCTION_ASSEMBLY",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_CATION_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_MICROTUBULE_POLYMERIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_PHOSPHOLIPID_METABOLIC_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_MRNA_PROCESSING",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CELL_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CELL_PROJECTION_ASSEMBLY",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_REGULATION_OF_MICROTUBULE_BASED_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_HEPATICOBILIARY_SYSTEM_DEVELOPMENT",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_REGULATION_OF_CYTOSKELETON_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways !=   "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_MRNA_PROCESSING",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_ION_TRANSPORT", ]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CELLULAR_COMPONENT_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_NEGATIVE_REGULATION_OF_NEURON_DEATH",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_MULTICELLULAR_ORGANISMAL_SIGNALING",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_TRANSPORT_ALONG_MICROTUBULE", ]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_DNA_REPLICATION",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_DNA_REPLICATION", ]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_INTRACELLULAR_PROTEIN_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_POSITIVE_REGULATION_OF_MITOTIC_CELL_CYCLE",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE_PHASE_TRANSITION",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_NEGATIVE_REGULATION_OF_PROTEIN_CONTAINING_COMPLEX_ASSEMBLY",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_CELL_PROJECTION_ASSEMBLY",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_RNA_PROCESSING",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_SYNAPSE_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_DEVELOPMENTAL_MATURATION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_BIOLOGICAL_PROCESS_INVOLVED_IN_INTERACTION_WITH_HOST",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_BODY_FLUID_LEVELS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_MACROAUTOPHAGY",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_EXPORT_FROM_CELL",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_ANATOMICAL_STRUCTURE_SIZE",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_MONOSACCHARIDE_BIOSYNTHETIC_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_ORGANIC_HYDROXY_COMPOUND_BIOSYNTHETIC_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_LOCALIZATION_WITHIN_MEMBRANE",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CELL_KILLING",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_HEMOSTASIS",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_NEGATIVE_REGULATION_OF_PROTEIN_POLYMERIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CEREBRAL_CORTEX_DEVELOPMENT",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_REGULATION_OF_TELOMERE_MAINTENANCE_VIA_TELOMERE_LENGTHENING",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_NEGATIVE_REGULATION_OF_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_CATION_TRANSMEMBRANE_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_SPINDLE_LOCALIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_STEROL_BIOSYNTHETIC_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_ISOPRENOID_METABOLIC_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_CELLULAR_COMPONENT_SIZE",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_MIDBRAIN_DEVELOPMENT",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOBP_ENDOPLASMIC_RETICULUM_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_ALCOHOL_BIOSYNTHETIC_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_REGULATION_OF_BLOOD_CIRCULATION",]

#group timepoint 2 pathways together
tp2_GOBP <- GO_merged[GO_merged$timepoint == "timepoint 2",] 
tp2_GOBP <-tp2_GOBP[order(tp2_GOBP$`FDR.q.val`),] 

tp3_GOBP <- GO_merged[GO_merged$timepoint == "timepoint 3",] 
tp3_GOBP <-tp3_GOBP[order(tp3_GOBP$`FDR.q.val`),] 

tp2_GOBP$Pathways
tp3_GOBP$Pathways

GO_merged$Pathways <- factor(GO_merged$Pathways, levels = c( "GOBP_AXON_DEVELOPMENT",                                            
 "GOBP_NEURON_PROJECTION_GUIDANCE",                                  
 "GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION",       
 "GOBP_SUBSTANTIA_NIGRA_DEVELOPMENT",                                
 "GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE",                       
 "GOBP_REGULATION_OF_AXONOGENESIS",                                  
 "GOBP_CELL_JUNCTION_ORGANIZATION",                                  
 "GOBP_CELL_PART_MORPHOGENESIS",                                     
 "GOBP_POSITIVE_REGULATION_OF_CELL_PROJECTION_ORGANIZATION",         
 "GOBP_REGULATION_OF_MICROTUBULE_CYTOSKELETON_ORGANIZATION",         
 "GOBP_NEURAL_NUCLEUS_DEVELOPMENT",                                  
 "GOBP_POTASSIUM_ION_TRANSPORT",                                     
 "GOBP_PROTEIN_LOCALIZATION_TO_CELL_PERIPHERY",                      
 "GOBP_POSITIVE_REGULATION_OF_DEVELOPMENTAL_GROWTH",                 
 "GOBP_PHOSPHOLIPID_BIOSYNTHETIC_PROCESS",                           
 "GOBP_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION",             
 "GOBP_REGULATION_OF_MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION",
 "GOBP_VESICLE_MEDIATED_TRANSPORT_TO_THE_PLASMA_MEMBRANE",           
 "GOBP_NEURON_MIGRATION",                                            
 "GOBP_DEVELOPMENTAL_GROWTH_INVOLVED_IN_MORPHOGENESIS",              
 "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",                      
 "GOBP_CARBOHYDRATE_BIOSYNTHETIC_PROCESS",                           
 "GOBP_PROTEIN_LOCALIZATION_TO_PLASMA_MEMBRANE",                     
 "GOBP_REGULATION_OF_PROTEIN_POLYMERIZATION",                        
 "GOBP_NEUROTRANSMITTER_SECRETION",                                  
 "GOBP_ESTABLISHMENT_OF_CELL_POLARITY",                              
 "GOBP_NEURON_DEVELOPMENT",                                          
 "GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY",                 
 "GOBP_NEGATIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS",             
 "GOBP_SYNAPTIC_SIGNALING",                                          
 "GOBP_NEGATIVE_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION",
"GOBP_AXONAL_TRANSPORT",                                            
"GOBP_MACROAUTOPHAGY",                                              
"GOBP_CYTOSKELETON_DEPENDENT_INTRACELLULAR_TRANSPORT",              
"GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",                                 
"GOBP_NEURON_PROJECTION_EXTENSION",                                 
"GOBP_LIPID_BIOSYNTHETIC_PROCESS",  
"GOBP_ACTIN_FILAMENT_POLYMERIZATION",                               
"GOBP_AXO_DENDRITIC_TRANSPORT",                                  
"GOBP_MICROTUBULE_BASED_MOVEMENT",                                  
"GOBP_REGULATED_EXOCYTOSIS",                                        
"GOBP_ACTIN_POLYMERIZATION_OR_DEPOLYMERIZATION",                    
"GOBP_NEURON_DEATH",                                                
"GOBP_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY",                
"GOBP_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE",
 "GOBP_RECOMBINATIONAL_REPAIR",                                   
 "GOBP_DNA_DEPENDENT_DNA_REPLICATION",                            
 "GOBP_DNA_GEOMETRIC_CHANGE",                                     
 "GOBP_DNA_CONFORMATION_CHANGE",                                  
 "GOBP_RNA_SPLICING_VIA_TRANSESTERIFICATION_REACTIONS",           
 "GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE",                         
 "GOBP_DNA_RECOMBINATION",                                        
 "GOBP_MRNA_METABOLIC_PROCESS",                                   
 "GOBP_RNA_SPLICING",                                             
 "GOBP_NCRNA_PROCESSING",                                         
 "GOBP_DNA_METABOLIC_PROCESS",                                    
 "GOBP_DNA_REPAIR",                                               
 "GOBP_DOUBLE_STRAND_BREAK_REPAIR",                               
 "GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_BY_RNA_POLYMERASE_II",
 "GOBP_CHROMATIN_REMODELING"))

plot_NES_merged(GO_merged, GO_db ="GOBP", firsttimepoint = "2", secondtimepoint = "3")


#############do pathway analysis for DCX timepoint 1 vs DCX timepoint 2

###############for filtering GSEA results
#for DCX timepoint1: FDR.q.val` <= 0.25
#for DCX timepoint2: FDR.q.val` <= 0.25
#in total there were 427 pathways which is too many


#adjust to FDR.q.value cutoff DCX of 0.1 for tp 1
#FDR.q.value cutoff of DCX of 0.05 for DCX tp 2

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/")


#read in GSEA results
filter_gsea_DCX_tps(gsea_report_na_pos="results/dcx_tp2_vs_tp1_GO_human_all_pathways.Gsea.1657654730497/gsea_report_for_tp1_1657654730497.tsv", 
                    gsea_report_for_na_neg="results/dcx_tp2_vs_tp1_GO_human_all_pathways.Gsea.1657654730497/gsea_report_for_tp2_1657654730497.tsv",
                    timepoint1 = "timepoint 1", timepoint2="timepoint 2")

#GOCC
#filter for GOCC:
extract_specific_GO_DCX (filtered_DCX=filtered_DCX, filtered_DCX2=filtered_DCX2,
                         GO_db="GOCC", timepoint1="timepoint 1", timepoint2="timepoint 2")

#keep only FDR.qvalue < 0.042 for DCX timepoint 2; keep it more stringent to fit on plot
DCX2_GO <-DCX2_GO[DCX2_GO$FDR.q.val<= 0.042, ]

GO_merged <-rbind(DCX_GO, DCX2_GO)

GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_SPLICEOSOMAL_TRI_SNRNP_COMPLEX",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_ACTIN_FILAMENT_BUNDLE",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOCC_PLASMA_MEMBRANE_PROTEIN_COMPLEX",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_SIDE_OF_MEMBRANE",]
GO_merged <- GO_merged[GO_merged$Pathways !=  "GOCC_VACUOLAR_MEMBRANE",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_CELL_SURFACE",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_ACTIN_FILAMENT_BUNDLE",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_SIDE_OF_MEMBRANE",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_VACUOLAR_MEMBRANE",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_U2_SNRNP",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_U1_SNRNP",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_U12_TYPE_SPLICEOSOMAL_COMPLEX",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_U2_TYPE_CATALYTIC_STEP_2_SPLICEOSOME",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_SPLICEOSOMAL_COMPLEX",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOCC_VACUOLE",]

#group timepoint 1 pathways together
GOCC_tp1 <- GO_merged[GO_merged$timepoint == "timepoint 1",]
GOCC_tp1$Pathways

GOCC_tp2 <- GO_merged[GO_merged$timepoint == "timepoint 2",]
GOCC_tp2$Pathways

GO_merged$Pathways <- factor(GO_merged$Pathways, levels = c(
  "GOCC_NUCLEAR_PROTEIN_CONTAINING_COMPLEX", "GOCC_U2_TYPE_SPLICEOSOMAL_COMPLEX",       "GOCC_NUCLEAR_BODY",                      
  "GOCC_SPLICEOSOMAL_SNRNP_COMPLEX",         "GOCC_NUCLEAR_SPECK",                      "GOCC_SM_LIKE_PROTEIN_FAMILY_COMPLEX",    
  "GOCC_PRECATALYTIC_SPLICEOSOME",           "GOCC_CATALYTIC_STEP_2_SPLICEOSOME",       "GOCC_CAJAL_BODY",                        
  "GOCC_METHYLTRANSFERASE_COMPLEX",          "GOCC_CHROMATIN",                          "GOCC_RIBONUCLEOPROTEIN_COMPLEX",         
  "GOCC_CYTOPLASMIC_STRESS_GRANULE",         "GOCC_P_BODY",                             "GOCC_HISTONE_METHYLTRANSFERASE_COMPLEX", 
  "GOCC_CHROMOSOME",                         "GOCC_ATPASE_COMPLEX",                     "GOCC_TRANSFERASE_COMPLEX",               
  "GOCC_TRANSCRIPTION_REGULATOR_COMPLEX",    "GOCC_NUCLEOLUS",                          "GOCC_CHROMOSOME_TELOMERIC_REGION",
  "GOCC_INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE",                     "GOCC_PLASMA_MEMBRANE_REGION",                                    
  "GOCC_CELL_CELL_JUNCTION",                                         "GOCC_BASOLATERAL_PLASMA_MEMBRANE",                               
  "GOCC_BASAL_PART_OF_CELL",                                         "GOCC_EXTRINSIC_COMPONENT_OF_CYTOPLASMIC_SIDE_OF_PLASMA_MEMBRANE",
  "GOCC_AXON",                                                       "GOCC_LEADING_EDGE_MEMBRANE",                                     
  "GOCC_CELL_LEADING_EDGE",                                          "GOCC_MEMBRANE_MICRODOMAIN",                                      
  "GOCC_RECYCLING_ENDOSOME",                                         "GOCC_ACTIN_FILAMENT",                                            
  "GOCC_RUFFLE",                                                     "GOCC_CELL_PROJECTION_MEMBRANE",                                  
  "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",                   "GOCC_NEURON_PROJECTION",                                         
  "GOCC_ACTIN_CYTOSKELETON",                                         "GOCC_ADHERENS_JUNCTION",                                         
  "GOCC_APICAL_PART_OF_CELL",                                        "GOCC_SUPRAMOLECULAR_POLYMER",                                    
  "GOCC_PHAGOCYTIC_VESICLE",                                         "GOCC_ENDOCYTIC_VESICLE",                                         
  "GOCC_CELL_CORTEX",                                                "GOCC_POLYMERIC_CYTOSKELETAL_FIBER",                              
  "GOCC_APICAL_PLASMA_MEMBRANE",                                     "GOCC_ANCHORING_JUNCTION",                                        
  "GOCC_EXTERNAL_ENCAPSULATING_STRUCTURE",                           "GOCC_ACTIN_BASED_CELL_PROJECTION",                               
  "GOCC_MICROTUBULE_ASSOCIATED_COMPLEX",                             "GOCC_FILOPODIUM",                                                
  "GOCC_SYNAPSE",                                                    "GOCC_NEURON_SPINE",                                              
  "GOCC_POSTSYNAPSE",                                                "GOCC_TRANS_GOLGI_NETWORK",                                       
  "GOCC_ACTOMYOSIN",                                                 "GOCC_NEURON_TO_NEURON_SYNAPSE",                                  
  "GOCC_ENDOPLASMIC_RETICULUM_LUMEN",                                "GOCC_DENDRITIC_TREE",                                            
  "GOCC_CELL_SUBSTRATE_JUNCTION"))                                   

plot_NES_merged(GO_merged, GO_db ="GOCC", firsttimepoint = "1", secondtimepoint = "2")


#filter for GOMF:
extract_specific_GO_DCX (filtered_DCX=filtered_DCX, filtered_DCX2=filtered_DCX2,
                         GO_db="GOMF", timepoint1="timepoint 1", timepoint2="timepoint 2")

#group tp1 pathways together; group tp2 pathways together
DCX_GO$Pathways
DCX2_GO$Pathways

GO_merged$Pathways <- factor(GO_merged$Pathways, levels = c("GOMF_HELICASE_ACTIVITY",                                                  
 "GOMF_ATP_DEPENDENT_ACTIVITY_ACTING_ON_RNA",                               
 "GOMF_HISTONE_BINDING",                                                    
 "GOMF_CATALYTIC_ACTIVITY_ACTING_ON_A_NUCLEIC_ACID",                        
 "GOMF_TRANSCRIPTION_COREGULATOR_ACTIVITY",                                 
 "GOMF_CATALYTIC_ACTIVITY_ACTING_ON_RNA",                                   
 "GOMF_TRANSCRIPTION_REGULATOR_ACTIVITY",                                   
 "GOMF_NUCLEASE_ACTIVITY",                                                  
 "GOMF_CHROMATIN_BINDING",                                                  
 "GOMF_MRNA_BINDING",                                                       
 "GOMF_TRANSCRIPTION_COACTIVATOR_ACTIVITY",                                 
 "GOMF_TRANSCRIPTION_COREPRESSOR_ACTIVITY",                                 
 "GOMF_MODIFICATION_DEPENDENT_PROTEIN_BINDING",                             
"GOMF_MRNA_3_UTR_BINDING",                                                 
 "GOMF_PRE_MRNA_BINDING",                                                   
"GOMF_TRANSCRIPTION_FACTOR_BINDING",                                       
"GOMF_ATP_DEPENDENT_ACTIVITY_ACTING_ON_DNA",                               
 "GOMF_CATALYTIC_ACTIVITY_ACTING_ON_DNA",                                   
"GOMF_SINGLE_STRANDED_RNA_BINDING",                                        
"GOMF_ATP_HYDROLYSIS_ACTIVITY",                                            
"GOMF_SNRNA_BINDING",                                                      
"GOMF_DNA_BINDING_TRANSCRIPTION_FACTOR_BINDING",                           
"GOMF_NUCLEOSOME_BINDING",                                                 
"GOMF_NUCLEAR_RECEPTOR_BINDING",                                           
"GOMF_RIBONUCLEOPROTEIN_COMPLEX_BINDING",                                  
"GOMF_HYDROLASE_ACTIVITY_ACTING_ON_ESTER_BONDS",                           
 "GOMF_NUCLEOSOMAL_DNA_BINDING",                            
 "GOMF_PROTEIN_METHYLTRANSFERASE_ACTIVITY",                                 
 "GOMF_RIBONUCLEASE_ACTIVITY",                                              
 "GOMF_RNA_POLYMERASE_II_SPECIFIC_DNA_BINDING_TRANSCRIPTION_FACTOR_BINDING",
"GOMF_DNA_HELICASE_ACTIVITY",  
"GOMF_STRUCTURAL_CONSTITUENT_OF_CYTOSKELETON",                        "GOMF_GTPASE_ACTIVITY",                                              
 "GOMF_INTEGRIN_BINDING",                                              "GOMF_INORGANIC_MOLECULAR_ENTITY_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
 "GOMF_CATION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",                     "GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",                       
 "GOMF_LIPID_BINDING",                                                 "GOMF_ACTIVE_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",                
 "GOMF_GUANYL_NUCLEOTIDE_BINDING",                                     "GOMF_PHOSPHOLIPID_BINDING",                                         
 "GOMF_PASSIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",                    "GOMF_GDP_BINDING",                                                  
 "GOMF_CARBOXYLIC_ACID_BINDING",                                       "GOMF_CALMODULIN_BINDING",                                           
 "GOMF_ACTIN_BINDING",                                                 "GOMF_CYTOSKELETAL_PROTEIN_BINDING",                                 
 "GOMF_ORGANIC_ACID_BINDING"))

plot_NES_merged(GO_merged, GO_db ="GOMF", firsttimepoint = "1", secondtimepoint = "2")

#for DCX tp 1 vs DCX tp 2 GOBP
extract_specific_GO_DCX (filtered_DCX=filtered_DCX, filtered_DCX2=filtered_DCX2,
                         GO_db="GOBP", timepoint1="timepoint 1", timepoint2="timepoint 2")

#filter GOBP results for DCX timepoint 1 by FDR.q.value <0.05
DCX_GO<-DCX_GO[DCX_GO$FDR.q.val <= 0.05, ]
DCX2_GO<-DCX2_GO[DCX2_GO$FDR.q.val <= 0.03, ]

GO_merged <- rbind(DCX_GO, DCX2_GO)

#remove redundant pathways
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_MRNA_PROCESSING",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_RNA_PROCESSING",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_REGULATION_OF_CALCIUM_ION_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_REGULATION_OF_ION_TRANSMEMBRANE_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_CALCIUM_ION_TRANSMEMBRANE_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_METAL_ION_TRANSPORT",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_CELL_PROJECTION_ASSEMBLY",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_CELL_JUNCTION_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways !="	GOBP_CELL_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_STRIATED_MUSCLE_CELL_DIFFERENTIATION",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_CELL_PART_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_ACTIN_FILAMENT_BASED_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_REGULATION_OF_ACTIN_FILAMENT_BASED_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_REGULATION_OF_ACTIN_FILAMENT_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_ACTIN_FILAMENT_ORGANIZATION",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_CELL_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_REGULATION_OF_ANATOMICAL_STRUCTURE_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_ACTIN_FILAMENT_BASED_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_REGULATION_OF_ACTIN_FILAMENT_BASED_PROCESS",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_REGULATION_OF_ACTIN_FILAMENT_LENGTH",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_ANATOMICAL_STRUCTURE_FORMATION_INVOLVED_IN_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways != "GOBP_CELLULAR_COMPONENT_MORPHOGENESIS",]
GO_merged <- GO_merged[GO_merged$Pathways !="GOBP_CHROMATIN_ORGANIZATION",]


plot_NES_merged(GO_merged, GO_db ="GOBP", firsttimepoint = "1", secondtimepoint = "2")
