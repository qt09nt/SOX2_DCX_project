setwd("C:\\Users\\qt09n\\Desktop\\Technical Analyst I UHN May 4 2021\\organoid group\\Sofia\\late prenatal human neurodevelopment resolved by single nucleus transcriptomics\\supplemental")

#prenatal ages of samples
age <- c("17 weeks gestation",
         "17 weeks gestation",
         "18 weeks gestation",
         "18 weeks gestation",
         "19 weeks gestation",
         "19 weeks gestation",
         "19.7 weeks gestation",
         "19.7 weeks gestation",
         "20.3 weeks gestation",
         "20.3 weeks gestation",
         "21.3 weeks gestation",
         "21.3 weeks gestation",
         "22 weeks gestation",
         "22 weeks gestation",
         "23.3 weeks gestation", 
         "23.3 weeks gestation", 
         "23.7 weeks gestation",
         "23.7 weeks gestation",
         "24.7 weeks gestation",
         "24.7 weeks gestation",
         "25.8 weeks gestation",
         "25.8 weeks gestation",
         "26 weeks gestation",
         "26 weeks gestation",
         "32.8 weeks gestation",
         "32.8 weeks gestation",
         "33.2 weeks gestation",
         "33.2 weeks gestation",
         "38.3 weeks gestation",
         "38.3 weeks gestation")

prenatal_sex <- c( "female",
                   "female",
                   "male",
                   "male",
                   "male",
                   "male",
                   "female",
                   "female",
                   "female",
                   "female",
                   "male",
                   "male",
                   "male",
                   "male",
                   "male",
                   "male",
                   "male",
                   "male",
                   "male",
                   "male",
                   "female",
                   "female",
                   "female",
                   "female",
                   "female",
                   "female",
                   "male",
                   "male",
                   "female",
                   "female"
)


cortical_plate <- read.csv("GSE217511_CorticalPlate_Seuratmetadata.csv")

unique(cortical_plate$celltypes)
# [1] "L4/5a CPN" "L2/3 CPN"  "IN"        "nIPC"      "L6 CPN"    "AC"        "UD"        "SPN"      
# [9] "gIPC"      "OPC"       "TAC"       "CRN"       "MG"        "BVC"  
#14 different cell types


germinal_matrix <- read.csv("GSE217511_GerminalMatrix_Seuratmetadata.csv")

unique(germinal_matrix$celltypes)
#[1] "UD"        "IN"        "L2/3 CPN"  "gIPC"      "L4/5a CPN" "TAC"       "nIPC"      "AC"   
#[9] "RG/AC"     "MSN"       "OPC"       "SPN"       "L6 CPN"    "MG"        "BVC"       "CP" 

subcortical_plate<-read.csv("GSE217511_subCorticalPlate_Seuratmetadata.csv")

unique(subcortical_plate$celltypes)
#[1] "AC-p"   "OPC"    "gIPC-A" "AC-f"   "gIPC"   "nIPC"   "gIPC-O" "EN"     "IN"     "preOL" 

subgerminal_matrix <- read.csv("GSE217511_subGerminalMatrix_Seuratmetadata.csv")
unique(subgerminal_matrix$celltypes)

# [1] "RG/AC"  "nIPC"   "gIPC"   "AC"     "oRG"    "mIPC"   "OPC"    "IN"     "EPD"    "gIPC-O" "preOL" 
# [12] "gIPC-A" "EN"     "tRG"    "CP"   

prenatal_celltypes_all <-c("L4/5a CPN",
                           "L2/3 CPN",  "IN",        "nIPC",      "L6 CPN",    "AC",        "UD",        "SPN",      
"gIPC",      "OPC",       "TAC",       "CRN",       "MG",        "BVC",  
"UD",        "IN",        "L2/3 CPN",  "gIPC",      "L4/5a CPN", "TAC",       "nIPC",      "AC",   
"RG/AC",     "MSN",       "OPC",       "SPN",       "L6 CPN",    "MG",        "BVC",       "CP", 
"AC-p",   "OPC",    "gIPC-A", "AC-f",   "gIPC",   "nIPC",   "gIPC-O", "EN",     "IN",     "preOL", 
"RG/AC",  "nIPC",   "gIPC",   "AC",     "oRG",    "mIPC",   "OPC",    "IN",     "EPD",    "gIPC-O", "preOL", 
"gIPC-A", "EN",     "tRG",    "CP" 
)

unique(prenatal_celltypes_all)

adult_cortex <- read.csv("GSE217511_Cortex_Seuratmetadata.csv")
unique(adult_cortex$celltypes)
#[1] "OL"        "OPC"       "L2/3 CPN"  "L4 CPN"    "MG"        "AC-f"      "IN 5HTR3a" "L5/6 CPN"  "AC-p"      "IN SOM"   
#[11] "IN PV"     "BVC"       "UD"        "preOL"    

adult_SVZ_caudate <- read.csv("GSE217511_SVZ_Caudate_Seuratmetadata.csv")
unique(adult_SVZ_caudate$celltypes)
#[1] "OL"        "MSN"       "MG"        "OPC"       "AC-p"      "UD"        "IN PV"     "IN 5HTR3a" "AC-f"      "IN SOM"   
#[11] "T CELL"    "EN"        "BVC"       "EPD" 

adult_all <- c("OL",        "OPC",       "L2/3 CPN",  "L4 CPN",    "MG",        "AC-f",      "IN 5HTR3a", "L5/6 CPN",
               "AC-p",      "IN SOM",   
               "IN PV",     "BVC",       "UD",        "preOL",  
               "OL",        "MSN",       "MG",        "OPC",       "AC-p",      "UD",        "IN PV",   
               "IN 5HTR3a", "AC-f",      "IN SOM",   
               "T CELL",    "EN",        "BVC",       "EPD"  )

unique(adult_all)
# [1] "OL"        "OPC"       "L2/3 CPN"  "L4 CPN"    "MG"        "AC-f"      "IN 5HTR3a" "L5/6 CPN"  "AC-p"      "IN SOM"    "IN PV"     "BVC"      
# [13] "UD"        "preOL"     "MSN"       "T CELL"    "EN"        "EPD"     

