library(tibble)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(tidyr)

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Enrichr/turquoise (neural differentiation) module list")

#read in turquoise module list genes (neural differentiation )
bioplanet <- read.table("BioPlanet_2019_table.txt", header=T, sep ="\t")
#bioplanet <- read.table("BioPlanet_2019_table.txt", header=T, row.names=1, sep ="\t")
go_biological_process <- read.table("GO_Biological_Process_2021_table.txt", header=T, row.names=1, sep ="\t")
GO_Cellular_Component<- read.table("GO_Cellular_Component_2021_table.txt", header=T, row.names=1, sep ="\t")

GO_Molecular_Function <- read.table("GO_Molecular_Function_2021_table.txt", header = T, row.names=1, sep="\t")

KEGG <- read.table("KEGG_2021_Human_table.txt", header = T, row.names=1, sep="\t")
Reactome <- read.table("Reactome_2016_table.txt", header = T, row.names=1, sep="\t")
WikiPathway <-  read.table("WikiPathway_2021_Human_table.txt", header = T, row.names=1, sep="\t")

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

processing(WikiPathway)

#Reactome needs further pre-processing
#Bioplanet also needs further pre-processing

go_biological_process_p <- processed
go_cellular_component_p <- processed
go_molecular_function_p <- processed
KEGG_p <- processed
WikiPathway_p <- processed

#not done yet
bioplanet_p <- processed
reactome_p <- processed

############################### Bioplanet results pre-processing

#extract Bioplanet row 167 which has several rows jumbled up together in one row
bioplanet_167 <- bioplanet[167,]


#data is the dataframe containing the extracted weird row from original dataframe
#original is the name of the dataframe with the original results output from Enrichr
#extractedrow is the row number of the extracted row in the or
preprocessing <- function(data, original, extractedrow) {
  row.names(data)<- data[,1]
  rows <- strsplit(rownames(data), "\n") 
  print(rows)
  rows <- unlist(rows)
  rows <- as.data.frame(rows)
  row_all <- data.frame()
  numrows <-nrow(rows)
  
  print(numrows)
  row.names(rows)<-rows[,1]
  #split the row by the tab character and place output of that row into the "all" dataframe
  for (i in 1:numrows){
    row <- c(strsplit(row.names(rows)[i], "\t")) %>% unlist()
    row_all <- rbind(row_all, row)
  
  }

  View(row_all)
  row.names(row_all)<- row_all[,1]
  column_names <- colnames(original)
  row.names(original)<-original[,1]
  colnames(row_all) <- column_names
  row_all <<- row_all
  last <- row.names(row_all)[numrows]
  
  row.names(original)[extractedrow] <- last  #rename row in original dataframe with the last row name of row_all
  original <- rbind(original, row_all)
  lastnum<-nrow(original)
  realrows <- lastnum - 1 
  original <- original[1:realrows,]        #get rid of the last row which is just repetition of the lastrow row name
  original <- original[,2:9]
  original$P.value <- as.numeric(original$P.value)
  original <<- original    #original Enrichr results table dataframe
}

preprocessing(bioplanet_167, bioplanet, 167)

original_473 <- original[473,]
rows <- strsplit(rownames(original_473), "\n") 
rows <- unlist(rows)
rows<-as.data.frame(rows)

row_all <- data.frame()

for (i in 1:278){
  row <- c(strsplit(rows[i,1], "\t")) %>% unlist()
  row_all <- rbind(row_all, row)
  }

View(row_all)
row.names(row_all) <- row_all[,1]
row_all <- row_all[2:9]
colnames(row_all)<-colnames(Reactome)
last<- row.names(row_all)[3]
row.names(reactome2)[130]<-last  #rename original bioplanet dataframe as Parkinson's disease
bioplanet <- rbind(original, row_all)
bioplanet <- bioplanet[1:714,]
bioplanet$P.value <- as.numeric(bioplanet$P.value)

processing(bioplanet)

bioplanet_p <-processed

########################################## Pre-processing Reactome
reactome247 <- Reactome[247, ]
rows<- strsplit(rownames(reactome247), "\n") 
rows<- as.data.frame(unlist(rows))
reactome <- Reactome[1:246,]

#split up the rows dataframe by the tab delimiter and place in columns
row_all2 <- data.frame()
for (i in 1:278){
  row <- c(strsplit(row.names(rows)[i], "\t")) %>% unlist()
  row_all2 <- rbind(row_all2, row)
}

row.names(row_all2)<-row_all2[,1]
row_all2<- row_all2[,2:9]
colnames(row_all2) <-colnames(Reactome)
reactome <- rbind(reactome, row_all2)

################pre-process reactome row 129 as well

reactome129 <- reactome[129,]
rows<- strsplit(rownames(reactome129), "\n") 
rows<- as.data.frame(unlist(rows))

row_all6 <- data.frame()
row.names(rows)<- rows[,1]
row.names(row_all6)<-row_all6[,1]
for (i in 1:196){
  row <- c(strsplit(row.names(rows)[i], "\t")) %>% unlist()
  row_all6 <- rbind(row_all6, row)
}

rows<- row_all6[, 2:9]
colnames(rows)<-colnames(reactome)
last<- row.names(rows)[196]
row.names(reactome)[129]<-last
rows <- rows[1:195,]
reactome <- rbind(reactome, rows)


reactome130<- reactome[130,]
rows<- strsplit(rownames(reactome130), "\n") 
rows<- as.data.frame(unlist(rows))
row_all7 <- data.frame()

row.names(rows)<- rows[,1]

for (i in 1:3){
  row <- c(strsplit(row.names(rows)[i], "\t")) %>% unlist()
  row_all7 <- rbind(row_all7, row)
}
row.names(row_all7)<- row_all7[,1]
rows<- row_all7[, 2:9]
colnames(rows)<-colnames(reactome)
last <- row.names(rows)[3]
row.names(reactome)[130]<-last
rows<- rows[1:2,]
reactome<- rbind(reactome, rows)


reactome131<- reactome[131,]
rows<- strsplit(rownames(reactome131), "\n") 
rows<- as.data.frame(unlist(rows))
row_all8 <- data.frame()

for (i in 1:52){
  row <- c(strsplit(row.names(rows)[i], "\t")) %>% unlist()
  row_all8 <- rbind(row_all8, row)
}
row.names(row_all8)<- row_all8[,1]
row_all8<-row_all8[2:9]
colnames(row_all8)<-colnames(reactome)
last<-row.names(row_all8)[52]

row.names(reactome)[131]<-last
rows<-row_all8[1:51,]
reactome<-rbind(reactome, rows)

reactome$P.value<- as.numeric(reactome$P.value)
processing(reactome)

reactome_p <- processed
#########################
write_csv(WikiPathway_p, "WikiPathway_p.csv")
write_csv(bioplanet_p, "bioplanet_p.csv")
write_csv(go_biological_process_p, "go_biological_process_p.csv")
write_csv(go_cellular_component_p, "go_cellular_component_p.csv")
write_csv(go_molecular_function_p, "go_molecular_function_p.csv")
write_csv(KEGG_p, "KEGG_p.csv")
write_csv(reactome_p, "reactome_p.csv")
