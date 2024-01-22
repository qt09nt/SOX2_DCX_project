rm(list = ls())

# change this to the directory you are working in
setwd("~/R projects/2020-2021/2021-06-29_AIM1_Analysis_Queenie/")

setwd("C:/Users/queenie.tsang/Desktop/USP15/2021-06-29_AIM1_Analysis_Queenie/2021-06-29_AIM1_Analysis_Queenie")

# run code containing additional functions I made
source("perseus_functions.R")
source("main_analysis_functions.R")

# load data (already pre-processed in Perseus)
data = read.table("data/data_preprocessed_in_perseus.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)
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

write.table(data_pos, "results/data_pos.txt", sep = "\t", quote = F)
write.table(metadata_pos, "results/metadata_pos.txt", sep = "\t", quote = F, row.names = F)

### Ugi's data
data_ugi = read.table("data/Ugi/data_brain.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)
metadata_ugi = read.table("data/ugi/metadata_brain.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)
rownames(data_ugi) = get_gene_names(data_ugi, gene_colname = "Gene", protein_colname = "Protein")
data_ugi = select_columns(data_ugi, "1", "32")
colnames(data_ugi) = paste(colnames(data_ugi), metadata_ugi$`Sample Type`)

write.table(data_ugi, "results/data_ugi_processed.txt", sep = "\t", quote = F)

### Jen's data
data_jen = read.table("data/Jen/data.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)
metadata_jen = read.table("data/Jen/metadata.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)
rownames(data_jen) = get_gene_names(data_jen, gene_colname = "Gene", protein_colname = "Protein")
data_jen = select_columns(data_jen, "1", "48")

NP = metadata_jen$`Neurological Region` %in% c("DSVZ", "LGE", "MGE")
MN = metadata_jen$`Neurological Region` %in% c("CX", "STRIA")

metadata_jen$group[NP] = "NP"
metadata_jen$group[MN] = "MN"
metadata_jen$group[(!MN)&(!NP)] = "other"
metadata_jen_pooled = metadata_jen[NP|MN,]

data_jen_pooled = data_jen[, NP|MN]
colnames(data_jen_pooled) = paste(colnames(data_jen_pooled), metadata_jen_pooled$`Neurological Region`)

write.table(data_jen_pooled, "results/data_jen_pooled.txt", sep = "\t", quote = F)
write.table(metadata_jen_pooled, "results/metadata_jen_pooled.txt", sep = "\t", quote = F)

write.table(data_jen, "results/data_jen.txt", sep ="\t", quote = F,
            col.names = T, row.names = T)
write.table(metadata_jen, "results/metadata_jen.txt", sep ="\t", col.names = T,
            row.names = F, quote = F)

