library(tibble)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(tidyr)

## Input data comes from Enrichr, after inputing the list of red (transcription regulation) module genes

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Enrichr/red transcriptional module")

#read in red module list genes (transcription regulation )
bioplanet <- read.table("BioPlanet_2019_table.txt", header=T, row.names=1, sep ="\t")
go_biological_process <- read.table("GO_Biological_Process_2021_table.txt", header=T, row.names=1, sep ="\t")
GO_Cellular_Component<- read.table("GO_Cellular_Component_2021_table.txt", header=T, row.names=1, sep ="\t")

GO_Molecular_Function <- read.table("GO_Molecular_Function_2021_table.txt", header = T, row.names=1, sep="\t")

KEGG <- read.table("KEGG_2021_Human_table.txt", header = T, row.names=1, sep="\t")
Reactome <- read.table("Reactome_2016_table.txt", header = T, row.names=1, sep="\t")
WikiPathway <-  read.table("WikiPathway_2021_Human_table.txt", header = T, row.names=1, sep="\t")

#filter for significant results (p-value < 0.05)
bioplanet_sig <- subset(bioplanet, P.value < 0.05)

#get rid of trailing empty NA columns:
bioplanet_sig[,9:17]<- NULL

new <- bioplanet_sig

new <- new %>% mutate(term = rownames(new))

#split Genes column into separate rows
new <- new %>% 
  mutate(Genes = strsplit(as.character(Genes), ";")) %>% 
  unnest(Genes)

new_wider <- new %>% 
  pivot_wider(names_from = term, values_from = P.value)

#replace NA with 0
new_wider[is.na(new_wider)] <- 0

#keep only terms and gene info columns
bp <- new_wider[,7:79]

#group rows by gene names
bp<- bp %>%
  group_by(Genes) %>%
  summarise_each((sum))

bp[bp == 0] <- NA

bioplanet_processed <- bp

#For GO_Molecular_Function dataframe
#some weird thing going on with row 3 where the input file is not quite separating by tab delimiters properly so the other rows 
#are all being read into the row name of row 3 instead of across the columns
GO_Molecular_Function_row3 <- GO_Molecular_Function[3,] 

GO_Mol_Func <- strsplit(rownames(GO_Molecular_Function_row3), "\n")  #split by newline character
GO_Mol_Func <- unlist(GO_Mol_Func)

all <- data.frame()

#split the row by the tab character and place output of that row into the "all" dataframe
for (i in 1:8){  # 8 is the number of rows
  row <- c(strsplit(GO_Mol_Func[i], "\t")) %>% unlist()
  all <- rbind(all, row)
}

column_names <- colnames(GO_Molecular_Function)
row.names(all) <- all[,1]
all <- all[,2:9]
colnames(all) <- column_names
last <- rownames(all) 

#rename original GO Molecular Function 
row.names(GO_Molecular_Function)[3] <- last[8]
all <- all[1:7,]
GO_Molecular_Function <- rbind(GO_Molecular_Function, all)
GO_Molecular_Function[,2] <- as.numeric(GO_Molecular_Function[,2])

#######################################################################

##Pre-process Reactome dataframe, because row 7 is all jumbled together into 1 row name

reactome_row7 <- Reactome[7,]
reactome_row <- strsplit(rownames(reactome_row7), "\n")  #split by newline character
reactome_row <- unlist(reactome_row)

reactome_all <- data.frame()
#split the row by the tab character and place output of that row into the "all" dataframe
for (i in 1:21){
  row <- c(strsplit(reactome_row[i], "\t")) %>% unlist()
  reactome_all <- rbind(reactome_all, row)
}

column_names <- colnames(Reactome)
row.names(reactome_all) <- reactome_all[,1]
reactome_all <- reactome_all[,2:9]
colnames(reactome_all) <- column_names
last <- rownames(reactome_all) 

#rename original Reactome dataframe
row.names(Reactome)[7]<- last[21]
reactome_all <- reactome_all[1:20,]
reactome <- rbind(Reactome, reactome_all)
reactome[,2]<- as.numeric(reactome[,2])

#####################pre-process GO Biological Processes data frame because row 140 consists of several rows jumbled into it
#keep all columns which are not NA
go_biological_process <- go_biological_process[,1:8]

#extract GO Biological Processes row 140:
go_biological_process_row140 <- go_biological_process[140,] 
go_bio_proc140 <- strsplit(rownames(go_biological_process_row140), "\n")  #split by newline character
go_bio_proc140 <- unlist(go_bio_proc140)

all_gbp <- data.frame()

#split the row by the tab character and place output of that row into the "all" dataframe
for (i in 1:473){  # 473 is the number of rows
  row <- c(strsplit(go_bio_proc140[i], "\t")) %>% unlist()
  all_gbp <- rbind(all_gbp, row)
}

row.names(all_gbp) <- all_gbp[,1]  #set column 1 as rownames 
all_gbp <- all_gbp[,2:9]
column_names <- colnames(go_biological_process)
colnames(all_gbp)<- column_names


go_biological_process <- go_biological_process[1:139,] #remove the last row of original GO Biological Processes dataframe 

#bind together go_biological_process and all_gbp dataframes:
go_biological_process <- rbind(go_biological_process, all_gbp)

go_biological_process$P.value <- as.numeric(go_biological_process$P.value)

#go_biological_process row 45 also has some problems where several rows are jammed into row 45 instead of separating by tab delimiter
#extract row 45:
gbp_row45 <- go_biological_process[45,]

go_bio_proc45 <- strsplit(rownames(gbp_row45), "\n")  #split by newline character
go_bio_proc45 <- unlist(go_bio_proc45)

all_gbp <- data.frame()

#split the row by the tab character and place output of that row into the "all" dataframe
for (i in 1:6){  # 6 is the number of rows
  row <- c(strsplit(go_bio_proc45[i], "\t")) %>% unlist()
  all_gbp <- rbind(all_gbp, row)
}

row.names(all_gbp) <- all_gbp[,1]  #set column 1 as rownames 
all_gbp <- all_gbp[,2:9]
column_names <- colnames(go_biological_process)
colnames(all_gbp)<- column_names
all_gbp$P.value <- as.numeric(all_gbp$P.value)

last<- rownames(all_gbp)[6]
row.names(go_biological_process)[45]<-last  #rename row 45 of original dataframe with proper row name

all_gbp <- all_gbp[1:5,]
go_biological_process <- rbind(go_biological_process, all_gbp)

###################################################################
#function for processing ENRICHR output dataframes 
processing <- function(data){
  data <- subset(data, P.value < 0.05)
  data[,9:17]<- NULL
  #data <- data[,1:8]  #keep only first 8 columns 
  data <- data %>% mutate(term = rownames(data))
  data <- data %>% 
    mutate(Genes = strsplit(as.character(Genes), ";")) %>% 
    unnest(Genes)
  new_wider <- data %>% 
    pivot_wider(names_from = term, values_from = P.value)
  new_wider[is.na(new_wider)] <- 0
  new_wider <- new_wider[, 7:ncol(new_wider)]
  new_wider <- new_wider %>%
    group_by(Genes) %>%
    summarise_each((sum))
  new_wider[new_wider == 0] <- NA
  processed <<- as.data.frame(new_wider)
}

#run processing function on input dataframe:
processing(go_biological_process)

#save processed Enrichr results
go_biological_process_p <- processed
go_cellular_component_p <- processed
GO_Molecular_Function_p <- processed
KEGG_p <- processed

#need to pre-process the Reactome dataframe before running the processing function:
Reactome_p <- processed
WikiPathway_p <- processed

#write to csv:
write.csv(bioplanet_processed, 'bioplanet_p.csv') 
write.csv(go_biological_process_p, 'go_biological_process_p.csv') 
write.csv(go_cellular_component_p, 'go_cellular_component_p.csv')
write.csv(GO_Molecular_Function_p, 'go_molecular_function_p.csv')
write.csv(KEGG_p, 'KEGG_p.csv')
write.csv(Reactome_p, 'reactome_p.csv')
write.csv(WikiPathway_p, 'WikiPathway_p.csv')

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Enrichr/red transcriptional module/processed")

bioplanet <- read.csv(file = "bioplanet_p.csv", row.names = 1)
go_biological_process <- read.csv(file = "go_biological_process_p.csv", row.names = 1)
go_cellular_component <- read.csv(file = "go_cellular_component_p.csv", row.names = 1)
go_molecular_function <- read.csv(file = "go_molecular_function_p.csv", row.names = 1)
KEGG <- read.csv(file = "KEGG_p.csv", row.names = 1)
reactome <- read.csv(file = "reactome_p.csv", row.names = 1)
wikipathway <- read.csv(file = "WikiPathway_p.csv", row.names = 1)
