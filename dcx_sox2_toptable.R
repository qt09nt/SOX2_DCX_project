#this code is for finding the genes which are upregulated in SOX2 and the ones upregulated in DCX
#based on the DCX vs SOX2 volcano plot cutoff

setwd("/Users/queenietsang/Desktop/UHN Technical Analyst/Kevin_Sangster_2021-06-29_AIM1_Analysis_Queenie/")

# run code containing additional functions I made
source("perseus_functions.R")
source("main_analysis_functions.R")

DCX_SOX2_top_table = read.table("results/DCX_vs_SOX2_top_table.txt", sep = "\t", quote = "")

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
  
top_table <- DCX_SOX2_top_table

top_table$differentially_expressed <- "NA"

top_table<-top_table %>% mutate(DE = case_when(
  logFC > 1.333 & adj.P.Val < 0.1 ~ "DCX",
  logFC < 1.333 & adj.P.Val < 0.1 ~ "SOX2"
))

top_table$differentially_expressed <- NULL

newdata <- mydata[ which(mydata$gender=='F' 
                         & mydata$age > 65), ]

SOX2 <- top_table[which(top_table$DE == "SOX2"), ]
DCX <- top_table[which(top_table$DE == "DCX"), ]  

write.table(SOX2, "SOX2.txt", sep= "\t")
write.table(DCX, "DCX.txt", sep="\t")
