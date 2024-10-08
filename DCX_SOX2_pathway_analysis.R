#pathway analysis comparing Sofia's all DCX+ to all SOX2+ 

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/")

# run code containing additional functions Kevin Sangster made
source("R/perseus_functions.R")
source("R/main_analysis_functions.R")
source("R/protein_organoid_working_matrix_main_analysis_functions.R")

library(tidyverse)
library(viridis)

library(forcats)    #used for handling factors for gsea pathways dataframe & reordering by FDR.qvalue

#load the data (working matrix)
data = read.table("data/WorkingmatrixR.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)
rownames(data) = get_gene_names(data)

# save protein to gene mapping separately
gene_protein_map = cbind.data.frame(protein_name = data$`Protein name`, gene = rownames(data))
gene_protein_map = gene_protein_map[gene_protein_map$protein_name !="",]
write.table(gene_protein_map, "results/gene_protein_map.txt", sep = "\t", col.names = T, row.names = F)

# remove unused columns
data = select_columns(data, "DCX_YFP(-)_tp1_s.1", "SOX2_YFP(+)_tp2_s.3")
data = as.matrix(data)

# create annotations for timepoint and cell type
timepoints = rep(NA, ncol(data))
timepoints[grepl("tp1", colnames(data))] = 1
timepoints[grepl("tp2", colnames(data))] = 2
timepoints[grepl("tp3", colnames(data))] = 3

cell_type = rep(NA, ncol(data))
cell_type[grepl("DCX", colnames(data))] = "DCX"
cell_type[grepl("SOX2", colnames(data))] = "SOX2"

YFP = rep(NA, ncol(data))
YFP[grepl("YFP\\(\\+\\)", colnames(data))] = "+"
YFP[grepl("YFP\\(\\-\\)", colnames(data))] = "-"

metadata = data.frame(sample = colnames(data), tp = timepoints, cell_type = cell_type, YFP = YFP, stringsAsFactors = F)

# remove SOX2 tp3 (outlier)
data = data[,!(metadata$cell_type == "SOX2" & metadata$tp == "tp3")]
metadata = metadata[!(metadata$cell_type == "SOX2" & metadata$tp == "tp3"),]

# create separate data frame containing only the YFP+
data_pos = data[,metadata$YFP == "+"]
metadata_pos = metadata[metadata$YFP == "+",]

# save everything
write.table(data, "results/data_processed.txt", sep = "\t", quote = F)
write.table(metadata, "results/metadata.txt", sep = "\t", quote = F, row.names = F)

write.table(data_pos, "results/data_pos.txt", sep = "\t", quote = F, row.names = T)
write.table(metadata_pos, "results/metadata_pos.txt", sep = "\t", quote = F, row.names = F)



#### for the pathway analysis Run GSEA was used

#the files used for GSEA

#data_pos.gct contains the working matrix data which has been filtered for YFP+ samples only, 
#excluding SOX2+ YFP+ timepoint 3 data

#groups.cls.txt is the file which contains the group info ie. whether the sample is DCX or SOX2
#Human_AllPathways_June01_2017_symbol.gmt was used for database search

#Note: for "Metric for ranking genes" parameter; the tTest selection kept giving "java.lang.NullPointerException error" 
#message so "Signal2Noise" was selected instead



#load in the GSEA results
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/results/DCX_vs_SOX2.Gsea.1642445010245")

#filter GSEA results comparing DCX to SOX2
#returns filtered results for DCX as "filtered_DCX" dataframe
#return filtered GSEA results for SOx2 as "filtered_SOX2" dataframe

#set the `NOM.p.val` cutoff to <= 0.05
#set the FDR.q.val` cutoff to <= 0.05

#when NOM p. val cutoff is 0.1, and FDR.q.val cutoff is 0.1, the filtered DCX dataframe still has
#222 pathways and SOX2 filtered dataframe still has 242 pathways which makes it too crowded to plot
filter_gsea(gsea_report_na_pos = 'gsea_report_for_DCX_1642445010245.tsv',
                     gsea_report_for_na_neg = 'gsea_report_for_SOX2_1642445010245.tsv', cell_type1 = "DCX", cell_type2 = "SOX2")


#plot GSEA NES scores for DCX
plot_NES(cell_type="DCX", filtered_df = filtered_DCX)

plot_NES(cell_type="SOX2", filtered_df = filtered_SOX2)


#trying to split databases into 2 groups for GSEA analysis:
#first group = Human All Pathways June 01_2017_symbols.gmt, Hallmark, KEGG

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/results/")

filter_gsea(gsea_report_na_pos = 'DCX_vs_SOX_human_allpathways_hallmark_KEGG.Gsea.1642465992178/gsea_report_for_DCX_1642465992178.tsv',
            gsea_report_for_na_neg = 'SOX2_vs_DCX_human_allpathways_hallmark_KEGG.Gsea.1642469335354/gsea_report_for_SOX2_1642469335354.tsv', cell_type1 = "DCX", cell_type2 = "SOX2")

#plot GSEA NES scores for DCX
plot_NES(cell_type="DCX", filtered_df = filtered_DCX)


#SOX2 enriched pathway results filtered by:
#`NOM.p.val` <= 0.1
#`FDR.q.val` <= 0.25 
plot_NES(cell_type="SOX2", filtered_df = filtered_SOX2)


########## include all databases:
filter_gsea(gsea_report_na_pos = 'DCX_vs_SOX2_all_datasets.Gsea.1642471512606/gsea_report_for_DCX_1642471512606.tsv',
            gsea_report_for_na_neg = 'SOX2_vs_DCX_all_datasets.Gsea.1642471791186/gsea_report_for_SOX2_1642471791186.tsv', cell_type1 = "DCX", cell_type2 = "SOX2")

#ERROR: discrete value entered when continuous value expected
plot_NES(cell_type="DCX", filtered_df = filtered_DCX)


##split where databases includes: reactome, GO_biological_processes, GO_molecular_functions, and
#GO_cellular_components, and Human_allPathways
filter_gsea(gsea_report_na_pos = 'DCX_vs_SOX2_reactome_goBP_GOMF_GOCC_human_allpathways.Gsea.1642472737217/gsea_report_for_DCX_1642472737217.tsv',
            gsea_report_for_na_neg = 'SOX2_vs_DCX_reactome_goBP_GOMF_GOCC_human_allpathways.Gsea.1642472990776/gsea_report_for_SOX2_1642472990776.tsv', cell_type1 = "DCX", cell_type2 = "SOX2")

plot_NES(cell_type="DCX", filtered_df = filtered_DCX)
plot_NES(cell_type="SOX2", filtered_df = filtered_SOX2)

###NOTE:: way too crowded in terms of pathways that make it past the filters


################ split databases: include only REactome & Human AllPathways
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/results/")
filter_gsea(gsea_report_na_pos = 'DCX_vs_SOX_reactome_human_allPathways.Gsea.1642474222462/gsea_report_for_DCX_1642474222462.tsv',
            gsea_report_for_na_neg = 'SOX2_vs_DCX_reactome_human_allPathways.Gsea.1642474323463/gsea_report_for_SOX2_1642474323463.tsv', cell_type1 = "DCX", cell_type2 = "SOX2")

#for DCX:
#NOM.p.val` <= 0.1
#FDR.q.val` <= 0.1
plot_NES(cell_type="DCX", filtered_df = filtered_DCX)

#for SOX2:
#`NOM.p.val` <= 0.01,]
#`FDR.q.val` <= 0.03, ]
plot_NES(cell_type="SOX2", filtered_df = filtered_SOX2)


############# split databases include: Gene_ONtology_Biological_Processes, 
#GO_Molecular_Functions, GO_Cellular_Components
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/results/")

filter_gsea(gsea_report_na_pos = 'DCX_SOX2_gsea_GOBP_GOMF_GOCC.Gsea.1642476474401/gsea_report_for_DCX_1642476474401.tsv',
            gsea_report_for_na_neg = 'SOX2_vs_DCX_gsea_GOBP_GOMF_GOCC.Gsea.1642476682432/gsea_report_for_SOX2_1642476682432.tsv', cell_type1 = "DCX", cell_type2 = "SOX2")


#for DCX:
#NOM.p.val` <= 0.05
#FDR.q.val` <= 0.01
plot_NES(cell_type="DCX", filtered_df = filtered_DCX)

#for SOX2:
#`NOM.p.val` <= 0.1,]
#`FDR.q.val` <= 0.1, ]
plot_NES(cell_type="SOX2", filtered_df = filtered_SOX2)

#save this filtered DCX table and filtered SOX2 gsea pathways table
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/results")
write.csv(filtered_DCX, "filtered_DCX.csv")
saveRDS(filtered_DCX, "filtered_DCX.rds")
write.csv(filtered_SOX2, "filtered_SOX2.csv")
saveRDS(filtered_SOX2, "filtered_SOX2.rds")


###############################  Comparison of Common and Unique DCX & SOX2 pathways
#databases include: Gene_ONtology_Biological_Processes, 
#GO_Molecular_Functions, GO_Cellular_Components
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/results/")

#filter out any pathways which have `NOM.p.val` <= 0.1, and FDR.q.val` <= 0.1:
filter_gsea(gsea_report_na_pos = 'DCX_SOX2_gsea_GOBP_GOMF_GOCC.Gsea.1642476474401/gsea_report_for_DCX_1642476474401.tsv',
            gsea_report_for_na_neg = 'SOX2_vs_DCX_gsea_GOBP_GOMF_GOCC.Gsea.1642476682432/gsea_report_for_SOX2_1642476682432.tsv', cell_type1 = "DCX", cell_type2 = "SOX2")

#find the pathways which are unique to DCX 
DCX_pathways_unique <-unique.comparisons(dataframe1 = filtered_DCX, dataframe2 = filtered_SOX2)


#find the pathways which are unique to SOX2
SOX2_pathways_unique <- unique.comparisons(dataframe1 = filtered_SOX2, dataframe2 = filtered_DCX)
#output: all filtered SOX2 pathways were unique to SOX2

#find the pathways which are common between DCX and SOX2
common_pathways <- common.comparisons(dataframe1 = filtered_DCX, dataframe2 = filtered_SOX2)
#output: 0 common pathways when comparing the filtered DCX and SOX2 gsea pathway results


####### COmparisons for unique & common pathways between DCX & SOX2
#try filtering by a less stringent NOM p.value & re-try pathway comparisons:
#filter out any pathways which have `NOM.p.val` <= 0.25, and FDR.q.val` <= 0.25:
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/results/")

filter_gsea(gsea_report_na_pos = 'DCX_SOX2_gsea_GOBP_GOMF_GOCC.Gsea.1642476474401/gsea_report_for_DCX_1642476474401.tsv',
            gsea_report_for_na_neg = 'SOX2_vs_DCX_gsea_GOBP_GOMF_GOCC.Gsea.1642476682432/gsea_report_for_SOX2_1642476682432.tsv', cell_type1 = "DCX", cell_type2 = "SOX2")

#find the pathways which are unique to DCX 
DCX_pathways_unique <-unique.comparisons(dataframe1 = filtered_DCX, dataframe2 = filtered_SOX2)

#find the pathways which are unique to SOX2
SOX2_pathways_unique <- unique.comparisons(dataframe1 = filtered_SOX2, dataframe2 = filtered_DCX)

#find the pathways which are common between DCX and SOX2
common_pathways <- common.comparisons(dataframe1 = filtered_DCX, dataframe2 = filtered_SOX2)


################################ merging pathways in filtered DCX gsea results and filtered SOX2 results

#make a copy of filtered DCX pathways and SOX2 pathways
DCX_merged <- filtered_DCX
SOX2_merged <- filtered_SOX2

#delete the GS.DETAILS column which don't have any info
DCX_merged[,3] <- NULL
SOX2_merged[,3] <- NULL
#############For SOX2 pathways: 

###merge similar pathways together:
#merge GOMF_Helicase Activity & GOMF_DNA_Helicase_Activity
#(remove GOMF_Helicase Activity)
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOMF_HELICASE_ACTIVITY",]

#merge spliceosome pathways/keep just 1:
spliceosome <- subset(SOX2_merged, grepl("SPLICEOSOME", SOX2_merged$Pathways))
spliceosome$Pathways

#keep GOCC_U2_TYPE_CATALYTIC_STEP_2_SPLICEOSOME as the representative pathway for spliceosome
spliceosome <- spliceosome[spliceosome$Pathways != "GOCC_U2_TYPE_CATALYTIC_STEP_2_SPLICEOSOME",  ]
spliceosome$Pathways

rownames(SOX2_merged)<-SOX2_merged$Pathways

#remove the redundant spliceosome pathways from SOX2 dataframe
SOX2_merged <- SOX2_merged[(rownames(SOX2_merged) != spliceosome$Pathways[1:(length(spliceosome$Pathways))]),]
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOBP_DNA_REPLICATION",]

#merge "SPLICEOSOMAL" related pathways:
spliceosomal <- subset(SOX2_merged, grepl("SPLICEOSOMAL", SOX2_merged$Pathways))

#keep these pathways as representative:
spliceosomal <- spliceosomal[spliceosomal$Pathways != "GOBP_SPLICEOSOMAL_COMPLEX_ASSEMBLY", ]
spliceosomal <- spliceosomal[spliceosomal$Pathways != "GOCC_SPLICEOSOMAL_SNRNP_COMPLEX", ]
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOCC_SPLICEOSOMAL_COMPLEX",]


#remove redundant spliceosomal pathways:
SOX2_merged <- SOX2_merged[(rownames(SOX2_merged) != spliceosomal$Pathways[1:(length(spliceosomal$Pathways))]),]

#merge redundant snRNP pathways (small nuclear ribonucleoproteins - RNA-protein complexes that combine with 
#unmodified pre-mRNA and various other proteins to form a spliceosome, a large RNA-protein molecular 
#complex upon which splicing of pre-mRNA occurs)
SNRNP <- subset(SOX2_merged, grepl("SNRNP", SOX2_merged$Pathways))

#keep this SNRNP pathway as representative:
SNRNP <- SNRNP[SNRNP$Pathways != "GOCC_SPLICEOSOMAL_SNRNP_COMPLEX",]

#remove these other SNRNP pathways (redundant)
SOX2_merged <- SOX2_merged[(rownames(SOX2_merged) != SNRNP$Pathways[1]),]
SOX2_merged <- SOX2_merged[(rownames(SOX2_merged) != SNRNP$Pathways[2]),]
SOX2_merged <- SOX2_merged[(rownames(SOX2_merged) != SNRNP$Pathways[3]),]

#remove some of redundant nuclear related pathways
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOCC_NUCLEAR_PORE",]
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOCC_NUCLEAR_SPECK",]
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOCC_NUCLEAR_BODY",]

SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOBP_POSITIVE_REGULATION_OF_RNA_SPLICING",]
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOBP_MRNA_METABOLIC_PROCESS",]

write.csv(SOX2_merged, "SOX2_merged_jan_18_2022_630pm.csv")

SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOBP_CHROMATIN_ORGANIZATION",]

#remove mRNA redundant pathways:
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOBP_MRNA_PROCESSING",]
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOBP_RNA_SPLICING",]
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOBP_RIBONUCLEOPROTEIN_COMPLEX_SUBUNIT_ORGANIZATION",]



#GOCC_SM_LIKE_PROTEIN_FAMILY_COMPLEX
#A protein complex containing members of the Like-Sm family of proteins, 
#which includes both the Sm proteins and the Lsm proteins, and which generally 
#form hexameric or heptameric ring structures which bind to RNA. While some of these ring 
#complexes may form independently of RNA, many only form in association with their target RNA. 
#In addition to Lsm-family proteins, many of these complexes contain additional protein members.
#Members of this family of complexes include the snRNPs which comprise the majority of the spliceosome.
#Others are involved in the 5' to 3' degradation pathways of mRNAs in the cytoplasm and of unspliced 
#transcripts in the nucleus, as well as other diverse roles. [GOC:bhm, GOC:krc, PMID:19121818, PMID:27627834]
#http://www.gsea-msigdb.org/gsea/msigdb/cards/GOCC_SM_LIKE_PROTEIN_FAMILY_COMPLEX

SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOCC_PROTEIN_DNA_COMPLEX",]

# Chromatin is a substance within a chromosome consisting of DNA and protein. ... The major proteins 
# in chromatin are histones, which help package the DNA in a compact form that fits in the cell nucleus.
# Changes in chromatin structure are associated with DNA replication and gene expression.

#remove viral gene expression pathway due to  relatively high FDR q value:
SOX2_merged <- SOX2_merged[SOX2_merged$Pathways != "GOBP_VIRAL_GENE_EXPRESSION",]

write.csv(SOX2_merged, "SOX2_merged_Jan_18_2021_9PM.csv")
saveRDS(SOX2_merged, "SOX2_merged_Jan_18_2021_9PM.rds")

############### merging DCX pathways (removing redundant ones)
#remove very similar morphogenesis related pathways
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_DEVELOPMENTAL_GROWTH_INVOLVED_IN_MORPHOGENESIS",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION",]

#removing microtubule redundant pathways
DCX_merged <- DCX_merged[DCX_merged$Pathways != c("GOBP_MICROTUBULE_BASED_MOVEMENT"),]
DCX_merged <- DCX_merged[DCX_merged$Pathways != c("GOBP_REGULATION_OF_MICROTUBULE_POLYMERIZATION"), ]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOCC_AXON_CYTOPLASM",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOCC_NEURON_PROJECTION",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOCC_NEURON_PROJECTION_CYTOPLASM",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_NEURON_PROJECTION_EXTENSION",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOCC_MICROTUBULE",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_MICROTUBULE_BASED_TRANSPORT",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_POSITIVE_REGULATION_OF_PROTEIN_POLYMERIZATION", ]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_CELL_PART_MORPHOGENESIS", ]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_REGULATION_OF_PROTEIN_POLYMERIZATION", ]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_POSITIVE_REGULATION_OF_CELL_PROJECTION_ORGANIZATION",]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_CELLULAR_COMPONENT_MORPHOGENESIS", ]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_REGULATION_OF_ANATOMICAL_STRUCTURE_SIZE", ]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOBP_PROTEIN_POLYMERIZATION", ]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOMF_TUBULIN_BINDING", ]
DCX_merged <- DCX_merged[DCX_merged$Pathways != "GOCC_SUPRAMOLECULAR_POLYMER", ]

write.csv(DCX_merged, "DCX_merged.csv")
saveRDS(DCX_merged, "DCX_merged.rds")


########### read in the filtered dataframes for DCX and SOX2 gsea pathways (filtered for FDR qvalue & NOM pvalue)


filtered_DCX <- readRDS("results/filtered_DCX.rds")
filtered_SOX2<- readRDS("results/filtered_SOX2.rds")

###### read in previous DCX_merged dataframe
DCX_merged <- read.csv("results/DCX_merged.csv")
SOX2_merged <- readRDS("results/SOX2_merged_Jan_18_2021_9PM.rds")


#get rid of extra column "X" which is a byproduct of importing csv file into R:
DCX_merged[,1]<-NULL

merge_SOX2_DCX <- rbind(SOX2_merged, DCX_merged)
merge_SOX2_DCX$Pathways <- as.factor(merge_SOX2_DCX$Pathways)

#plot_NES_merged(merge_SOX2_DCX)

#order so SOX2 appears on left and DCX on right 
merge_SOX2_DCX$cell_type <- factor(merge_SOX2_DCX$cell_type, levels = c("SOX2", "DCX"))

#order the dots by ascending FDR.q value (lower FDR q value is higher on plot)
merge_SOX2_DCX_reordered <- merge_SOX2_DCX %>% mutate(Pathways = fct_reorder(Pathways, desc(FDR.q.val)))


#order the pathways on the y-axis so that the SOX2 pathways are together, and DCX pathways are together
merge_SOX2_DCX_reordered$Pathways <- factor(merge_SOX2_DCX_reordered$Pathways, levels = c("GOCC_SITE_OF_POLARIZED_GROWTH", "GOCC_MICROTUBULE_ASSOCIATED_COMPLEX",
                                                        "GOCC_DISTAL_AXON",
                                                        "GOCC_AXON",
                                                        "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY",
                                                        "GOBP_NEURON_PROJECTION_GUIDANCE",
                                                        "GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION",
                                                        "GOBP_AXON_EXTENSION",
                                                        "GOBP_AXON_DEVELOPMENT", 
                                                        "GOBP_AXO_DENDRITIC_TRANSPORT",
                                                        "GOMF_STRUCTURAL_CONSTITUENT_OF_CYTOSKELETON", 
                                                        "GOBP_NEURON_DEVELOPMENT",
                                                        "GOBP_AXONAL_TRANSPORT",
                                                        "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
                                                        "GOMF_MICROTUBULE_BINDING", 
                                                        "GOCC_CYTOPLASMIC_REGION",
                                                        "GOCC_POLYMERIC_CYTOSKELETAL_FIBER",
                                                        "GOBP_REGULATION_OF_EXTENT_OF_CELL_GROWTH",
                                                        "GOBP_CYTOSKELETON_DEPENDENT_INTRACELLULAR_TRANSPORT",
                                                        "GOBP_REGULATION_OF_MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION",
                                                        "GOBP_REGULATION_OF_AXONOGENESIS",
                                                        "GOCC_PLASMA_MEMBRANE_REGION",
                                                        "GOBP_IMPORT_INTO_CELL", 
                                                        "GOBP_CELL_PROJECTION_ORGANIZATION",
                                                        "GOBP_TRANSPORT_ALONG_MICROTUBULE", 
                                                        "GOBP_NEURON_DIFFERENTIATION",
                                                        "GOBP_REGULATION_OF_CELLULAR_COMPONENT_SIZE",
                                                        "GOBP_NEUROGENESIS",
                                                        "GOBP_ORGANIC_ACID_TRANSPORT", 
                                                        "GOBP_CELL_JUNCTION_ORGANIZATION",
                                                        "GOBP_TAXIS", 
                                                        "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT",
                                                        "GOBP_LOCOMOTION", 
                                                        "GOBP_CELL_RECOGNITION", 
                                                        "GOBP_NEGATIVE_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION",
                                                        "GOCC_CELL_BODY",  
                                                        "GOBP_REGULATION_OF_CELL_SIZE",  
                                                        "GOCC_CILIUM",
                                                        "GOBP_ACTIN_FILAMENT_BUNDLE_ORGANIZATION",
                                                        "GOCC_ACTIN_BASED_CELL_PROJECTION",
                                                        "GOBP_REGULATION_OF_AUTOPHAGY", 
                                                        "GOBP_MEMBRANE_FUSION", 
                                                        "GOBP_REGULATION_OF_CELL_SUBSTRATE_JUNCTION_ORGANIZATION",
                                                        "GOCC_FILOPODIUM",
                                                        "GOBP_NEURON_MIGRATION",
                                                        "GOBP_SPLICEOSOMAL_COMPLEX_ASSEMBLY",
                                                        "GOCC_SM_LIKE_PROTEIN_FAMILY_COMPLEX", 
                                                        "GOCC_SPLICEOSOMAL_SNRNP_COMPLEX",
                                                        "GOCC_NUCLEAR_PROTEIN_CONTAINING_COMPLEX",
                                                        "GOCC_U2_TYPE_CATALYTIC_STEP_2_SPLICEOSOME",
                                                        "GOBP_RNA_SPLICING_VIA_TRANSESTERIFICATION_REACTIONS",
                                                        "GOCC_NUCLEAR_CHROMOSOME",
                                                        "GOBP_MRNA_SPLICE_SITE_SELECTION",
                                                        "GOBP_DNA_DEPENDENT_DNA_REPLICATION",
                                                        "GOBP_CHROMATIN_REMODELING",
                                                        "GOBP_PROTEIN_DNA_COMPLEX_SUBUNIT_ORGANIZATION", 
                                                        "GOCC_MICROBODY",
                                                        "GOBP_NUCLEOSOME_ORGANIZATION",
                                                        "GOCC_ATPASE_COMPLEX", 
                                                        "GOBP_NUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHETIC_PROCESS",
                                                        "GOBP_MODULATION_BY_HOST_OF_SYMBIONT_PROCESS",
                                                        "GOBP_DNA_REPAIR",
                                                        "GOMF_ATP_DEPENDENT_ACTIVITY_ACTING_ON_DNA",
                                                        "GOMF_DNA_HELICASE_ACTIVITY", 
                                                        "GOMF_CATALYTIC_ACTIVITY_ACTING_ON_DNA",
                                                        "GOBP_RESPONSE_TO_NUTRIENT", 
                                                        "GOBP_BIOLOGICAL_PROCESS_INVOLVED_IN_INTERACTION_WITH_SYMBIONT",
                                                        "GOBP_IMPORT_INTO_NUCLEUS", 
                                                        "GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS",
                                                        "GOBP_FATTY_ACID_CATABOLIC_PROCESS", 
                                                        "GOBP_GLYCOSYL_COMPOUND_METABOLIC_PROCESS",
                                                        "GOBP_LIPID_MODIFICATION",
                                                        "GOBP_LIPID_OXIDATION",
                                                        "GOBP_CELLULAR_AMINO_ACID_CATABOLIC_PROCESS", 
                                                        "GOMF_NADP_BINDING",
                                                        "GOBP_CHROMOSOME_ORGANIZATION"))


ggplot(merge_SOX2_DCX_reordered, aes(x = cell_type, y = Pathways, size = NES, color = FDR.q.val )) +
  geom_point(alpha = 0.7) + #scatterplot
  scale_size(range = c(2, 5), name = "NES") +  #change the size of the points and create a name for the size legend
  scale_color_viridis(discrete = FALSE) 


################# try making separate plots with GO Molecular Functions, GO Biological Processes and GO 

library("stringr")
library("dplyr")

#GOBP
#subset for GOBP in DCX filtered dataframe

DCX_GOBP <- filter(filtered_DCX, grepl("GOBP", filtered_DCX$Pathways))
SOX2_GOBP <- filter(filtered_SOX2, grepl("GOBP", filtered_SOX2$Pathways))

GOBP_merged <- rbind(DCX_GOBP, SOX2_GOBP)
#order so SOX2 appears on left and DCX on right 
GOBP_merged$cell_type <- factor(GOBP_merged$cell_type, levels = c("SOX2", "DCX"))

#order the dots by ascending FDR.q value (lower FDR q value is higher on plot)
GOBP_merged <- GOBP_merged %>% mutate(Pathways = fct_reorder(Pathways, desc(FDR.q.val)))


GOBP_merged$Pathways <- factor(GOBP_merged$Pathways, levels = c("GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY",
  "GOBP_NEURON_PROJECTION_GUIDANCE",
  "GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION",
  "GOBP_AXON_EXTENSION",
  "GOBP_AXON_DEVELOPMENT",
  "GOBP_AXO_DENDRITIC_TRANSPORT",
  "GOBP_CELLULAR_COMPONENT_MORPHOGENESIS",
  "GOBP_NEURON_DEVELOPMENT",
  "GOBP_AXONAL_TRANSPORT",
  "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
  "GOBP_REGULATION_OF_PROTEIN_POLYMERIZATION",
  "GOBP_CELL_PART_MORPHOGENESIS",
  "GOBP_REGULATION_OF_EXTENT_OF_CELL_GROWTH",
  "GOBP_CYTOSKELETON_DEPENDENT_INTRACELLULAR_TRANSPORT",
  "GOBP_POSITIVE_REGULATION_OF_CELL_PROJECTION_ORGANIZATION",
  "GOBP_REGULATION_OF_MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION",
  "GOBP_NEURON_PROJECTION_EXTENSION",
  "GOBP_REGULATION_OF_AXONOGENESIS",
  "GOBP_IMPORT_INTO_CELL",
  "GOBP_PROTEIN_POLYMERIZATION",
  "GOBP_CELL_PROJECTION_ORGANIZATION",
  "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION",
  "GOBP_REGULATION_OF_MICROTUBULE_POLYMERIZATION",
  "GOBP_TRANSPORT_ALONG_MICROTUBULE",
  "GOBP_NEURON_DIFFERENTIATION",
  "GOBP_REGULATION_OF_CELLULAR_COMPONENT_SIZE",
  "GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION",
  "GOBP_NEUROGENESIS",
  "GOBP_MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION",
  "GOBP_ORGANIC_ACID_TRANSPORT",
  "GOBP_MICROTUBULE_BASED_TRANSPORT",
  "GOBP_DEVELOPMENTAL_GROWTH_INVOLVED_IN_MORPHOGENESIS",
  "GOBP_CELL_JUNCTION_ORGANIZATION",
  "GOBP_TAXIS",
  "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT",
  "GOBP_LOCOMOTION",
  "GOBP_CELL_RECOGNITION",
  "GOBP_NEGATIVE_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION",
  "GOBP_REGULATION_OF_CELL_SIZE",
  "GOBP_REGULATION_OF_ANATOMICAL_STRUCTURE_SIZE",
  "GOBP_ACTIN_FILAMENT_BUNDLE_ORGANIZATION",
  "GOBP_REGULATION_OF_AUTOPHAGY",
  "GOBP_MEMBRANE_FUSION",
  "GOBP_REGULATION_OF_CELL_SUBSTRATE_JUNCTION_ORGANIZATION",
  "GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT",
  "GOBP_MICROTUBULE_BASED_MOVEMENT",
  "GOBP_NEURON_MIGRATION",
  "GOBP_POSITIVE_REGULATION_OF_PROTEIN_POLYMERIZATION",
  "GOBP_SPLICEOSOMAL_COMPLEX_ASSEMBLY",
  "GOBP_RNA_SPLICING_VIA_TRANSESTERIFICATION_REACTIONS",
  "GOBP_MRNA_SPLICE_SITE_SELECTION",
  "GOBP_DNA_REPLICATION",
  "GOBP_DNA_DEPENDENT_DNA_REPLICATION",
  "GOBP_MRNA_PROCESSING",
  "GOBP_RNA_SPLICING",
  "GOBP_CHROMATIN_REMODELING",
  "GOBP_PROTEIN_DNA_COMPLEX_SUBUNIT_ORGANIZATION",
  "GOBP_CHROMATIN_ORGANIZATION",
  "GOBP_NUCLEOSOME_ORGANIZATION",
  "GOBP_NUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHETIC_PROCESS",
  "GOBP_MODULATION_BY_HOST_OF_SYMBIONT_PROCESS",
  "GOBP_DNA_REPAIR",
  "GOBP_MRNA_CIS_SPLICING_VIA_SPLICEOSOME",
  "GOBP_RESPONSE_TO_NUTRIENT",
  "GOBP_BIOLOGICAL_PROCESS_INVOLVED_IN_INTERACTION_WITH_SYMBIONT",
  "GOBP_IMPORT_INTO_NUCLEUS",
  "GOBP_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS",
  "GOBP_POSITIVE_REGULATION_OF_RNA_SPLICING",
  "GOBP_MRNA_METABOLIC_PROCESS",
  "GOBP_RIBONUCLEOPROTEIN_COMPLEX_SUBUNIT_ORGANIZATION",
  "GOBP_FATTY_ACID_CATABOLIC_PROCESS",
  "GOBP_GLYCOSYL_COMPOUND_METABOLIC_PROCESS",
  "GOBP_LIPID_MODIFICATION",
  "GOBP_LIPID_OXIDATION",
  "GOBP_SPLICEOSOMAL_SNRNP_ASSEMBLY",
  "GOBP_CELLULAR_AMINO_ACID_CATABOLIC_PROCESS",
  "GOBP_CHROMOSOME_ORGANIZATION",
  "GOBP_VIRAL_GENE_EXPRESSION"))

ggplot(GOBP_merged, aes(x = cell_type, y = Pathways, size = NES, color = FDR.q.val )) +
  geom_point(alpha = 0.7) + #scatterplot
  scale_size(range = c(2, 5), name = "NES") +  #change the size of the points and create a name for the size legend
  scale_color_viridis(discrete = FALSE) 














