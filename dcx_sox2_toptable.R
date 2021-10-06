#this code is for finding the genes which are upregulated in SOX2 and the ones upregulated in DCX
#based on the DCX vs SOX2 volcano plot cutoffs

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Kevin Sangster code/2021-06-29_AIM1_Analysis_Queenie/2021-06-29_AIM1_Analysis_Queenie/")

library(dplyr)

# run code containing additional functions I made
source("perseus_functions.R")
source("main_analysis_functions.R")

DCX_SOX2_top_table = read.table("results/DCX_vs_SOX2_top_table.txt", sep = "\t", header = TRUE)

data <- read.table("data/working_matrix.txt", sep = "\t", header= T, quote = "")

#DCX was set as the control when getting the top table

#function used for plotting DCX vs SOX2 volcano plot
# volcano_plot = function(top_table, title){
#   v = EnhancedVolcano(top_table, x="logFC", y= "adj.P.Val", lab=row.names(top_table), title = title,
#                       pCutoff = 0.1, FCcutoff = 1.333, 
#                       ylim = c(0, max(-log10(top_table$adj.P.Val), na.rm=TRUE)),
#                       col = c("grey30", "#A4D49C", "#4E80A5", "#F68D91"),
#                       colAlpha = 1.0)
#   pdf(paste("figures/", title, "volcano.pdf"), width = 8, height = 8)
#   print(v)

#make a copy of DCX_SOX2_top_table  
top_table <- DCX_SOX2_top_table

#create a new column in top_table called DE (differentially expressed) and fill in with values of 
#DCX or SOX2 depending on the cutoff values used for the volcano plots
top_table<-top_table %>% mutate(DE = case_when(
  logFC > 1.333 & adj.P.Val < 0.1 ~ "DCX",
  logFC < 1.333 & adj.P.Val < 0.1 ~ "SOX2"
))


#subset the rows of the top table for SOX2 and DCX in their own separate tables
SOX2 <- top_table[which(top_table$DE == "SOX2"), ]
DCX <- top_table[which(top_table$DE == "DCX"), ]  

write.table(SOX2, "SOX2.txt", sep= "\t")
write.table(DCX, "DCX.txt", sep="\t")

#get list of SOX2 genes
SOX2_genes <-row.names(SOX2)
class(SOX2_genes)

#get list of DCX genes
DCX_genes <- row.names(DCX)

#subset working matrix for SOX2 expressed genes
SOX2_DE <- subset(data, Gene.names %in% SOX2_genes)

#subset working matrix for DCX expressed genes
DCX_DE <- subset(data, Gene.names %in% DCX_genes)

#include up to AR on working matrix which is Unique Peptides

#remove AP ('Peptides') and AQ ('Razor + unique peptides') from matrix
SOX2_DE$Peptides <- NULL
DCX_DE$Peptides <- NULL
SOX2_DE$`Razor...unique.peptides` <- NULL
DCX_DE$`Razor...unique.peptides` <- NULL

#keeping just columns from AI(GOBP name) to AR:
SOX2_AI_AR <- SOX2_DE[,c(35:53)]
DCX_AI_AR <- DCX_DE[,c(35:53)]

#remove columns AS(Sequence coverage [%])  to AY (MS/MS count)
SOX2_AI_AR[,c(9:15)] <- NULL
DCX_AI_AR[,c(9:15)] <- NULL

#remove protein.IDs and remove Majority.protein.IDs
SOX2_AI_AR[,c(9,10)]<- NULL
DCX_AI_AR[,c(9,10)]<- NULL

write.table(SOX2_AI_AR, "results/SOX2_DE.txt",  row.names=FALSE, col.names=TRUE, quote=FALSE,  sep="\t")
write.table(DCX_AI_AR, "results/DCX_DE.txt",  row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

write.csv(SOX2_AI_AR, "results/SOX2_DE.csv", row.names= FALSE)
write.csv(DCX_AI_AR, "results/DCX_DE.csv", row.names= FALSE)


##############################

#re-do the p-values plot for red transcriptional module genes using the DCX_SOX2 top table 
#FDR < 0.1 (adj.P.Val)

column_names <- colnames(DCX_SOX2_top_table)
colnames(DCX_SOX2_top_table)<- c("genes", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/gene module lists")

#read in red (transciption regulation) module list genes
red <- read.table("red (transcriptional regulation) module list.txt") %>% unlist()
red <- as.data.frame(red)

#extract the rows of the DCX_SOX2_top_table which contain the genes in the red transcription module
red_only <- DCX_SOX2_top_table[DCX_SOX2_top_table$genes %in% red$red, ] 

#filter for the significantly differentially expressed proteins between SOX2 and DCX
#FDR adj.P.val < 0.1
significant <- red_only[red_only$`adj.P.Val` <= 0.1,]

library(ggplot2)

#re-plot adjusted p-values for red transcriptional module genes 
sp <- ggplot(data = red_only, aes(x=reorder(genes, adj.P.Val, mean, .desc=FALSE), y=adj.P.Val)) + #re-order x-axis by y-axis p-values mean
  geom_point() +
  xlab("gene") +
  ylab("Adjusted p-value") +
  theme_bw() + #set the background to white
  theme(panel.grid.major.x = element_blank(),  #removing the grey grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90)) +
  labs(title= "Proteins Data Pre-processed in Perseus",
       subtitle="Differentially Expressed genes from red transcriptional regulation module list - SOX2_YFP(+) vs. DCX_YFP(+)")
sp + geom_hline(yintercept=0.1, linetype="dashed") 
