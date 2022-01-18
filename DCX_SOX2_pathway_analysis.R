#pathway analysis comparing Sofia's all DCX+ to all SOX2+ 

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/")

# run code containing additional functions Kevin Sangster made
source("R/perseus_functions.R")
source("R/main_analysis_functions.R")
source("R/protein_organoid_working_matrix_main_analysis_functions.R")

library(tidyverse)
library(viridis)

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


################ split databases: include only REactome & Human ALlPathways
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




