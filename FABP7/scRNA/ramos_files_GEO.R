#automate downloading the 108 sample files from Ramos et al (2022) single cell RNA study


library(RCurl)
library(downloader)

url = "http://www.ncbi.nlm.nih.gov/geo/download/?acc="
setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/late prenatal human neurodevelopment resolved by single nucleus transcriptomics")

samples <- read.table(file = "samples.txt", sep = '\t', header = FALSE, fill = TRUE)

samples<-as.data.frame(samples)
samples[,1]<-substr(samples[,1], 1, 10)

samples[,2]<-c("_11C_", "_11G_", "_56C_", "_56G_","_34C_", "_34G_", "_30C_", "_30G_", "_63C_", "_63G_", "_23C_",
               "_23G_", "_62C_", "_62G_", 
               "_64C_",
               "_64G_", 
               "_3C_", 
               "_3G_", 
               "_60C_", 
               "_60G_",
               
               "_10C_", 
               
               "_10G_", 
               
               "_6C_", 
               "_6G_", 
               "_31C_",
               "_31G_", 
               "_2C_", 
               "_2G_", 
               "_20C_", "_20G_")


require(downloader)

#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6720862&format=file&file=GSM6720862%5F34C%5Fbarcodes%2Etsv%2Egz

for (i in 1:nrow(samples)) {
  filename = samples[i,1]
  id = samples[i, 2]

  download(url = (paste(url, filename, "&format=file&file=", filename, id, "barcodes.tsv.gz", sep = "")), destfile = (paste(getwd(), "/", filename, id, "_barcodes.tsv.gz", sep = "")), mode = "wb")
  download(url = (paste(url, filename, "&format=file&file=", filename, id,"features.tsv.gz", sep = "")),  destfile = (paste(getwd(), "/", filename, id, ".tsv.gz", sep = "")), mode = "wb")
  download(url = (paste(url, filename, "&format=file&file=", filename, id, "matrix.mtx.gz", sep = "")), destfile = (paste(getwd(), "/", filename, id, ".mtx.gz", sep = "")), mode = "wb")
  
}

