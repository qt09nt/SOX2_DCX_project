library(tibble)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(tidyr)

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Enrichr/blue (DNA replication) module list")

#read in blue module list genes (transcription regulation )
bioplanet <- read.table("BioPlanet_2019_table.txt", header=T, row.names=1, sep ="\t")
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

processing(Reactome)

WikiPathway_p <- processed
bioplanet_p <- processed
go_biological_process_p <- processed
go_cellular_component_p <- processed
go_molecular_function_p <- processed
KEGG_p <- processed
reactome_p <- processed

write_csv(WikiPathway_p, "WikiPathway_p.csv")
write_csv(bioplanet_p, "bioplanet_p.csv")
write_csv(go_biological_process_p, "go_biological_process_p.csv")
write_csv(go_cellular_component_p, "go_cellular_component_p.csv")
write_csv(go_molecular_function_p, "go_molecular_function_p.csv")
write_csv(KEGG_p, "KEGG_p.csv")
write_csv(reactome_p, "reactome_p.csv")
