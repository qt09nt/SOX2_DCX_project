#Analysis of control organoid proteins compared to schizophrenia organoid proteins

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/")

source("R/perseus_functions.R")
source("R/main_analysis_functions.R")
source("R/protein_organoid_working_matrix_main_analysis_functions.R")

library(VennDiagram)
library(tidyverse)
library(viridis)
library(forcats)    #used for handling factors for gsea pathways dataframe & reordering by FDR.qvalue
library(dplyr)
library(stringr)  #for subsetting gsea results for separate Gene Ontology Molecular Functions (GOMF), Gene Ontology Cellular Components (GOCC)
#and GOBP (Biological Processes) plots
library(pheatmap)
library(MCMCglmm)
library(ggplot2)
library(EnhancedVolcano)
library(limma)
library(factoextra)
library(dendextend)
library(purrr)
library(gridExtra)
library(cowplot)
library(hrbrthemes)
library(viridis)
library(WGCNA)
library(RColorBrewer)
library(ggsignif)
library(colorspace)
library(ggpubr)
library(cluster)
library(factoextra)
library(DEGreport)

#load GSEA results

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/schizophrenia/pathway analysis/organoid_proteins_ctr_vs_scz.Gsea.1672859077431")

filter_gsea(gsea_report_na_pos = 'gsea_report_for_CTR_1672859077431.tsv',
            gsea_report_for_na_neg = 'gsea_report_for_SCZ_1672859077431.tsv', cell_type1 = "CTR", cell_type2 = "SCZ", cell_type1_pvalue = 0.4, cell_type2_pvalue = 0.3)

#rename filtered_DCX to filtered_CTR and rename filtered_SOX2 to filtered_SCZ

filtered_CTR <- filtered_DCX
filtered_SCZ <- filtered_SOX2 

rm(filtered_DCX)
rm(filtered_SOX2)

#find the unique pathways for control organoids and unique pathways for schizophrenia organoids
control_unique <-unique.comparisons(dataframe1=filtered_CTR, dataframe2 = filtered_SCZ)
scz_unique <- unique.comparisons(dataframe1 = filtered_SCZ, dataframe2 = filtered_CTR)

#filter for GOBP only 
GOBP <-extract_specific_GO_db(filtered_DCX = filtered_CTR, filtered_SOX2=filtered_SCZ, "GOBP", celltype1= "CTR", celltype2="SCZ")

#count frequency of words in pathways 
#https://countwordsfree.com/

#remove some of the redundant pathways
GOBP<- GOBP[GOBP$Pathways !="GOBP_LAMELLIPODIUM_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_LAMELLIPODIUM_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_CARBOHYDRATE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CARBOHYDRATE_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CARBOHYDRATE_BIOSYNTHETIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_PROTEIN_BINDING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_CARBOHYDRATE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CARBOHYDRATE_DERIVATIVE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_POSITIVE_REGULATION_OF_CYTOSKELETON_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_ACTIN_FILAMENT_POLYMERIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_ACTIN_FILAMENT_DEPOLYMERIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_ACTIN_FILAMENT_BASED_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ACTIN_FILAMENT_DEPOLYMERIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ACTIN_FILAMENT_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ACTIN_FILAMENT_BASED_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_ACTIN_FILAMENT_LENGTH",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ORGANIC_HYDROXY_COMPOUND_BIOSYNTHETIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PYRIMIDINE_CONTAINING_COMPOUND_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_SULFUR_COMPOUND_BIOSYNTHETIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_CYTOSKELETON_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CORTICAL_CYTOSKELETON_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PROTEIN_LOCALIZATION_TO_CYTOSKELETON",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MUSCLE_CELL_MIGRATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MUSCLE_SYSTEM_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MUSCLE_CONTRACTION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_LAMELLIPODIUM_ASSEMBLY",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_DENDRITE_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_VESICLE_MEDIATED_TRANSPORT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_INTERMEDIATE_FILAMENT_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RESPONSE_TO_ALCOHOL",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_ALTERNATIVE_MRNA_SPLICING_VIA_SPLICEOSOME",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NCRNA_PROCESSING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RNA_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_RNA_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RECEPTOR_MEDIATED_ENDOCYTOSIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MITOTIC_CELL_CYCLE",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE_PHASE_TRANSITION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RESPONSE_TO_TOXIC_SUBSTANCE",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PROTEIN_KINASE_B_SIGNALING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_PROTEOLYSIS_INVOLVED_IN_PROTEIN_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RAS_PROTEIN_SIGNAL_TRANSDUCTION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MICROTUBULE_NUCLEATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CIRCULATORY_SYSTEM_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELLULAR_ALDEHYDE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELL_PROJECTION_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ORGANOPHOSPHATE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELL_SUBSTRATE_JUNCTION_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NUCLEOTIDE_PHOSPHORYLATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_PROTEIN_MODIFICATION_BY_SMALL_PROTEIN_CONJUGATION_OR_REMOVAL",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_OSTEOBLAST_DIFFERENTIATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_POSITIVE_REGULATION_OF_PROTEOLYSIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_PROTEOLYSIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_PROTEOLYSIS",]
GOBP<- GOBP[GOBP$Pathways !=  "GOBP_REGULATION_OF_PROTEIN_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MACROMOLECULE_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_POSITIVE_REGULATION_OF_PROTEIN_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_SMALL_MOLECULE_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ACTIN_POLYMERIZATION_OR_DEPOLYMERIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_LAMELLIPODIUM_ASSEMBLY",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_RNA_SPLICING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RNA_PROCESSING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_MRNA_PROCESSING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RNA_SPLICING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MRNA_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MRNA_PROCESSING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PROTEIN_MODIFICATION_BY_SMALL_PROTEIN_CONJUGATION_OR_REMOVAL",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PROTEIN_POLYMERIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_PROTEIN_POLYMERIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_HEART_MORPHOGENESIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_DENDRITE_MORPHOGENESIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_GERM_CELL_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_EPITHELIAL_CELL_DIFFERENTIATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MONOCARBOXYLIC_ACID_BIOSYNTHETIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PROTEOLYSIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_POSITIVE_REGULATION_OF_PROTEOLYSIS_INVOLVED_IN_PROTEIN_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELLULAR_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PROTEOLYSIS_INVOLVED_IN_PROTEIN_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_SUPRAMOLECULAR_FIBER_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MONOSACCHARIDE_BIOSYNTHETIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MONOSACCHARIDE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_SKIN_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_KERATINOCYTE_DIFFERENTIATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_MUSCLE_CELL_DIFFERENTIATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_DENDRITE_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_MITOTIC_CELL_CYCLE_PHASE_TRANSITION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_DETOXIFICATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_DENDRITIC_SPINE_MORPHOGENESIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ENDOSOMAL_TRANSPORT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ENDOCYTOSIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RETROGRADE_TRANSPORT_ENDOSOME_TO_GOLGI",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELLULAR_RESPONSE_TO_HYDROGEN_PEROXIDE",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_CELLULAR_COMPONENT_SIZE",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PURINE_NUCLEOSIDE_MONOPHOSPHATE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_CYTOSKELETON_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_HEPATICOBILIARY_SYSTEM_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CARDIAC_CHAMBER_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_LYTIC_VACUOLE_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_BLOOD_PRESSURE",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RIBOSOMAL_LARGE_SUBUNIT_BIOGENESIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_VACUOLE_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_HEART_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CARDIAC_MUSCLE_CELL_ACTION_POTENTIAL",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_PEPTIDASE_ACTIVITY",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELLULAR_AMIDE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_RNA_SPLICING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_VASCULAR_PROCESS_IN_CIRCULATORY_SYSTEM",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELLULAR_OXIDANT_DETOXIFICATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELL_CELL_JUNCTION_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_POSTSYNAPSE_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways !=  "GOBP_STEROL_BIOSYNTHETIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_POSITIVE_REGULATION_OF_PROTEIN_CONTAINING_COMPLEX_ASSEMBLY",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ANATOMICAL_STRUCTURE_HOMEOSTASIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_VASCULATURE_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PEPTIDE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CARBOHYDRATE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ESTABLISHMENT_OF_CELL_POLARITY",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ANATOMICAL_STRUCTURE_MATURATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_BINDING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_ORGANELLE_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RIBONUCLEOSIDE_MONOPHOSPHATE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELLULAR_PROCESS_INVOLVED_IN_REPRODUCTION_IN_MULTICELLULAR_ORGANISM",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CYTOKINESIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MITOTIC_CYTOKINESIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RIBONUCLEOSIDE_DIPHOSPHATE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_GLUCOSE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_CELL_PROJECTION_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_UROGENITAL_SYSTEM_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_RESPIRATORY_SYSTEM_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_PROTEIN_CONTAINING_COMPLEX_ASSEMBLY",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_SEX_DIFFERENTIATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PYRIDINE_CONTAINING_COMPOUND_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PURINE_CONTAINING_COMPOUND_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PROTEIN_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_POSITIVE_REGULATION_OF_SIGNALING",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_PROTEIN_STABILITY",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_TELENCEPHALON_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_MIDBRAIN_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_DEVELOPMENTAL_MATURATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_PROTEIN_LOCALIZATION_TO_NUCLEUS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_POSITIVE_REGULATION_OF_VIRAL_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ICOSANOID_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_INTRACELLULAR_SIGNAL_TRANSDUCTION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CYTOSOLIC_TRANSPORT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_INTRACELLULAR_TRANSPORT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_IMPORT_INTO_NUCLEUS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ORGANIC_HYDROXY_COMPOUND_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NEGATIVE_REGULATION_OF_PEPTIDASE_ACTIVITY",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_CELLULAR_COMPONENT_BIOGENESIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_INTRACELLULAR_SIGNAL_TRANSDUCTION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_NUCLEOBASE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELLULAR_RESPONSE_TO_TOXIC_SUBSTANCE",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_GLUTATHIONE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_FIBROBLAST_PROLIFERATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_REGULATION_OF_VIRAL_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PHAGOCYTOSIS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_CELLULAR_RESPONSE_TO_XENOBIOTIC_STIMULUS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_PROTEIN_MODIFICATION_BY_SMALL_PROTEIN_CONJUGATION",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_SMALL_MOLECULE_BIOSYNTHETIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways != "GOBP_ISOPRENOID_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways !=  "GOBP_POSITIVE_REGULATION_OF_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_REGULATION_OF_SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_CELLULAR_RESPONSE_TO_EXTRACELLULAR_STIMULUS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_POSITIVE_REGULATION_OF_ORGANELLE_ASSEMBLY",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_CELL_MATURATION",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_RESPONSE_TO_CYTOKINE",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_SUBSTRATE_ADHESION_DEPENDENT_CELL_SPREADING",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_ALCOHOL_BIOSYNTHETIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_RETINA_HOMEOSTASIS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_ORGANONITROGEN_COMPOUND_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_POSITIVE_REGULATION_OF_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_NUCLEOBASE_CONTAINING_SMALL_MOLECULE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_CELLULAR_KETONE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_REGULATION_OF_PROTEIN_POLYMERIZATION",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_NUCLEOSIDE_DIPHOSPHATE_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_METENCEPHALON_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_MONOSACCHARIDE_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_ENDOSOME_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_RESPONSE_TO_HYDROGEN_PEROXIDE",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_SULFUR_COMPOUND_METABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_REGULATION_OF_ANATOMICAL_STRUCTURE_SIZE",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_NEUROMUSCULAR_PROCESS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_NEURAL_NUCLEUS_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_EPITHELIAL_CELL_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_REGULATION_OF_SUPRAMOLECULAR_FIBER_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_LOCOMOTORY_BEHAVIOR",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_CYTOSKELETON_DEPENDENT_CYTOKINESIS",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_PROTEIN_POLYUBIQUITINATION",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_POSITIVE_REGULATION_OF_LAMELLIPODIUM_ORGANIZATION",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_SUBSTANTIA_NIGRA_DEVELOPMENT",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_ERBB_SIGNALING_PATHWAY",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_CELL_CELL_JUNCTION_ASSEMBLY",]
GOBP<- GOBP[GOBP$Pathways !="GOBP_NUCLEOSIDE_MONOPHOSPHATE_METABOLIC_PROCESS",]

#group all control pathways together and group all schizophrenia pathways together
CTR_only<- GOBP[GOBP$cell_type=="CTR",]
SCZ_only<- GOBP[GOBP$cell_type=="SCZ",]

CTR_only$Pathways
SCZ_only$Pathways


GOBP$Pathways <- factor(GOBP$Pathways, levels = c(
"GOBP_CELLULAR_MODIFIED_AMINO_ACID_METABOLIC_PROCESS",              
"GOBP_PRESYNAPTIC_ENDOCYTOSIS",                                     
"GOBP_REGULATION_OF_ACTIN_FILAMENT_ORGANIZATION",                   
"GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",                         
"GOBP_MONOCARBOXYLIC_ACID_METABOLIC_PROCESS",                       
"GOBP_DENDRITIC_SPINE_DEVELOPMENT",                                 
"GOBP_ESTABLISHMENT_OR_MAINTENANCE_OF_CELL_POLARITY",               
"GOBP_HINDBRAIN_DEVELOPMENT",                                       
"GOBP_ADP_METABOLIC_PROCESS",                                       
"GOBP_REGENERATION",                                                
"GOBP_PYRUVATE_METABOLIC_PROCESS",                                  
"GOBP_REGULATION_OF_NEURON_DIFFERENTIATION",                        
"GOBP_REGULATION_OF_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION",
"GOBP_POSITIVE_REGULATION_OF_PROTEIN_BINDING",                      
"GOBP_POSTSYNAPSE_ORGANIZATION",                                    
"GOBP_RESPONSE_TO_GROWTH_FACTOR",                                   
"GOBP_UNSATURATED_FATTY_ACID_METABOLIC_PROCESS",                    
"GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION",               
"GOBP_FATTY_ACID_BIOSYNTHETIC_PROCESS",                             
"GOBP_CEREBRAL_CORTEX_DEVELOPMENT",                                 
"GOBP_NADH_METABOLIC_PROCESS",                                      
"GOBP_ORGANIC_ACID_METABOLIC_PROCESS",                              
"GOBP_SYNAPTIC_VESICLE_RECYCLING",                                  
"GOBP_NEURON_PROJECTION_ORGANIZATION",                              
"GOBP_FATTY_ACID_METABOLIC_PROCESS",                                
"GOBP_MEMBRANE_DOCKING",                                            
"GOBP_AXON_EXTENSION",     
"GOBP_CYTOPLASMIC_TRANSLATION",                                                
"GOBP_RIBOSOME_BIOGENESIS",                                                    
"GOBP_RRNA_METABOLIC_PROCESS",                                                 
"GOBP_NCRNA_METABOLIC_PROCESS",                                                
"GOBP_PEPTIDE_BIOSYNTHETIC_PROCESS",                                           
"GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",                                   
"GOBP_NEGATIVE_REGULATION_OF_MRNA_METABOLIC_PROCESS",                          
"GOBP_AMIDE_BIOSYNTHETIC_PROCESS",                                             
"GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS",                            
"GOBP_RNA_SPLICING_VIA_TRANSESTERIFICATION_REACTIONS",                         
"GOBP_ATP_BIOSYNTHETIC_PROCESS",                                               
"GOBP_PROTON_TRANSMEMBRANE_TRANSPORT",                                         
"GOBP_NEGATIVE_REGULATION_OF_MRNA_PROCESSING",                                 
"GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",                               
"GOBP_REGULATION_OF_MRNA_METABOLIC_PROCESS",                                   
"GOBP_RNA_STABILIZATION",                                                      
"GOBP_OXIDATIVE_PHOSPHORYLATION",                                              
"GOBP_SODIUM_ION_TRANSMEMBRANE_TRANSPORT",                                     
"GOBP_RIBOSOME_ASSEMBLY",                                                      
"GOBP_REGULATION_OF_MRNA_SPLICING_VIA_SPLICEOSOME",                            
"GOBP_NEGATIVE_REGULATION_OF_RNA_METABOLIC_PROCESS",                           
"GOBP_DEVELOPMENT_OF_PRIMARY_SEXUAL_CHARACTERISTICS",                          
"GOBP_PROTEIN_DNA_COMPLEX_SUBUNIT_ORGANIZATION",                               
"GOBP_NEGATIVE_REGULATION_OF_NUCLEOBASE_CONTAINING_COMPOUND_METABOLIC_PROCESS",
"GOBP_ENDOPLASMIC_RETICULUM_ORGANIZATION",                                     
"GOBP_CHROMATIN_REMODELING",                                                   
"GOBP_MALE_SEX_DIFFERENTIATION",                                               
"GOBP_NEGATIVE_REGULATION_OF_GENE_EXPRESSION",                              
"GOBP_RNA_MODIFICATION",                                                       
 "GOBP_U2_TYPE_PRESPLICEOSOME_ASSEMBLY"))      

plot_NES(filtered_df = GOBP, plottitle = "Ramos et al (2022) dataset - Control vs Schizophrenia Organoids")  


###### filter for Gene Ontology Cell Components (GOCC) enriched pathways

#filter for GOBP only 
GOCC <-extract_specific_GO_db(filtered_DCX = filtered_CTR, filtered_SOX2=filtered_SCZ, "GOCC", celltype1= "CTR", celltype2="SCZ")

GOCC_CTR <- GOCC[GOCC$cell_type == "CTR", ]
GOCC_SCZ <- GOCC[GOCC$cell_type == "SCZ",]

#filter out redundant pathways to make plot more legible
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways !="GOCC_CLATHRIN_COATED_VESICLE_MEMBRANE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways !="GOCC_SIDE_OF_MEMBRANE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways !="GOCC_CYTOPLASMIC_SIDE_OF_MEMBRANE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways !="GOCC_EXTRINSIC_COMPONENT_OF_CYTOPLASMIC_SIDE_OF_PLASMA_MEMBRANE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways !="GOCC_EXTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways !="GOCC_PLASMA_MEMBRANE_REGION",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways !="GOCC_INTERMEDIATE_FILAMENT_CYTOSKELETON",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways !="GOCC_ACTIN_FILAMENT",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways !="GOCC_CORTICAL_ACTIN_CYTOSKELETON",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways != "GOCC_FICOLIN_1_RICH_GRANULE_LUMEN",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways != "GOCC_VESICLE_MEMBRANE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways != "GOCC_TRANSPORT_VESICLE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways != "GOCC_COATED_VESICLE_MEMBRANE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways != "GOCC_ENDOSOME_MEMBRANE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways != "GOCC_EXTERNAL_ENCAPSULATING_STRUCTURE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways != "GOCC_SECRETORY_GRANULE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways != "GOCC_VACUOLE",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways != "GOCC_AUTOPHAGOSOME",]
GOCC_CTR<- GOCC_CTR[GOCC_CTR$Pathways != "GOCC_PHAGOCYTIC_VESICLE",]

GOCC_SCZ<- GOCC_SCZ[GOCC_SCZ$Pathways != "GOCC_POLYSOMAL_RIBOSOME",]
GOCC_SCZ<- GOCC_SCZ[GOCC_SCZ$Pathways != "GOCC_SMALL_RIBOSOMAL_SUBUNIT",]
GOCC_SCZ<- GOCC_SCZ[GOCC_SCZ$Pathways != "GOCC_CYTOSOLIC_LARGE_RIBOSOMAL_SUBUNIT",]
GOCC_SCZ<- GOCC_SCZ[GOCC_SCZ$Pathways != "GOCC_RIBOSOME",]
GOCC_SCZ<- GOCC_SCZ[GOCC_SCZ$Pathways != "GOCC_CYTOSOLIC_SMALL_RIBOSOMAL_SUBUNIT",]
GOCC_SCZ<- GOCC_SCZ[GOCC_SCZ$Pathways != "GOCC_INTRINSIC_COMPONENT_OF_ORGANELLE_MEMBRANE",]
GOCC_SCZ<- GOCC_SCZ[GOCC_SCZ$Pathways != "GOCC_ENVELOPE",]
GOCC_SCZ<- GOCC_SCZ[GOCC_SCZ$Pathways != "GOCC_INTRINSIC_COMPONENT_OF_ORGANELLE_MEMBRANE",]


GOCC_combined<- rbind(GOCC_CTR, GOCC_SCZ)

#order pathways so that CTR pathways are together and SCZ pathways are together

GOCC_combined$Pathways <- factor(GOCC_combined$Pathways, levels = c(
    
  "GOCC_ACTIN_CYTOSKELETON",                      
  "GOCC_PODOSOME",                      
  "GOCC_ENDOPEPTIDASE_COMPLEX",        
  "GOCC_PROTEASOME_COMPLEX",   
  "GOCC_ACTOMYOSIN",                    
  "GOCC_CELL_CORTEX",                             
  "GOCC_ACTIN_FILAMENT_BUNDLE",      
  "GOCC_LAMELLIPODIUM",                   
  "GOCC_FICOLIN_1_RICH_GRANULE",                  
  "GOCC_CONTRACTILE_FIBER",           
  "GOCC_CORTICAL_CYTOSKELETON",   
  "GOCC_SUPRAMOLECULAR_POLYMER",                  
  "GOCC_MEMBRANE_COAT",        
  "GOCC_UBIQUITIN_LIGASE_COMPLEX",  
  "GOCC_CELL_LEADING_EDGE",                       
  "GOCC_VESICLE_LUMEN",       
  "GOCC_PEPTIDASE_COMPLEX",        
  "GOCC_TRANSPORT_VESICLE_MEMBRANE",    
  "GOCC_VESICLE_COAT",             
  "GOCC_GABA_ERGIC_SYNAPSE",            
  "GOCC_INTERCELLULAR_BRIDGE",                    
  "GOCC_POSTSYNAPTIC_MEMBRANE",      
  "GOCC_CYTOPLASMIC_SIDE_OF_PLASMA_MEMBRANE",     
  "GOCC_CLATHRIN_COATED_PIT",                     
  "GOCC_SITE_OF_DOUBLE_STRAND_BREAK",  
  "GOCC_POLYMERIC_CYTOSKELETAL_FIBER",      
  "GOCC_CELL_DIVISION_SITE",                      
  "GOCC_INTERMEDIATE_FILAMENT",        
  "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",
  "GOCC_NUCLEAR_PORE",                            
  "GOCC_LATE_ENDOSOME",                 
  "GOCC_RUFFLE_MEMBRANE",                    
  "GOCC_PROTEASOME_ACCESSORY_COMPLEX",            
  "GOCC_SYNAPTIC_MEMBRANE",             
  "GOCC_MYELIN_SHEATH",  
  "GOCC_CYTOSOLIC_RIBOSOME",                                    
  "GOCC_RIBOSOMAL_SUBUNIT",                                    
  "GOCC_LARGE_RIBOSOMAL_SUBUNIT",                               
  "GOCC_RIBONUCLEOPROTEIN_COMPLEX",                            
  "GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX",              
  "GOCC_ORGANELLE_INNER_MEMBRANE",                             
  "GOCC_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX",          
  "GOCC_NUCLEOLUS",                                            
  "GOCC_RESPIRASOME",                                           
  "GOCC_MITOCHONDRIAL_ENVELOPE",                               
  "GOCC_CATALYTIC_STEP_2_SPLICEOSOME",                          
  "GOCC_SPLICEOSOMAL_COMPLEX",                                 
  "GOCC_POLYSOME",                                              
  "GOCC_NUCLEAR_SPECK",                                        
  "GOCC_CAJAL_BODY",                                            
  "GOCC_INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE",               
  "GOCC_INTRINSIC_COMPONENT_OF_ENDOPLASMIC_RETICULUM_MEMBRANE", 
  "GOCC_NUCLEAR_BODY",                                         
  "GOCC_TRANSPORTER_COMPLEX",                                   
  "GOCC_U2_TYPE_SPLICEOSOMAL_COMPLEX",                         
  "GOCC_U2_SNRNP",                                              
  "GOCC_PRECATALYTIC_SPLICEOSOME",                             
  "GOCC_PRESYNAPTIC_ACTIVE_ZONE",                               
  "GOCC_CHROMOSOME"))

plot_NES(filtered_df = GOCC_combined, plottitle = "Ramos et al (2022) dataset - Control vs Schizophrenia Organoids - Gene Ontology Cell Components")  


#filter for Gene Ontology Molecular Functions (GoMF) pathways
GOMF <-extract_specific_GO_db(filtered_DCX = filtered_CTR, filtered_SOX2=filtered_SCZ, "GOMF", celltype1= "CTR", celltype2="SCZ")

GOMF_CTR <- GOMF[GOMF$cell_type=="CTR",]
GOMF_SCZ <- GOMF[GOMF$cell_type=="SCZ",]



GOMF_CTR<- GOMF_CTR[GOMF_CTR$Pathways != "GOMF_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_THE_ALDEHYDE_OR_OXO_GROUP_OF_DONORS_NAD_OR_NADP_AS_ACCEPTOR",]



GOMF_combined <- rbind(GOMF_CTR, GOMF_SCZ)

GOMF_combined$Pathways<- factor(GOMF_combined$Pathways, levels =c(
"GOMF_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_THE_ALDEHYDE_OR_OXO_GROUP_OF_DONORS",                        
"GOMF_GTPASE_ACTIVITY",                                                                              
"GOMF_PROTEASE_BINDING",                                                                             
"GOMF_GUANYL_NUCLEOTIDE_BINDING",                                                                    
"GOMF_ACTIN_BINDING",                                                                                
"GOMF_CARBOXYLIC_ACID_BINDING",                                                                      
"GOMF_STRUCTURAL_CONSTITUENT_OF_CYTOSKELETON",                                                       
"GOMF_PROTEIN_HOMODIMERIZATION_ACTIVITY",                                                            
"GOMF_PEPTIDASE_REGULATOR_ACTIVITY",                                                                 
"GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT",                                                  
"GOMF_PHOSPHOLIPID_BINDING",                                                                         
"GOMF_LIPID_BINDING",                                                                                
"GOMF_ACTIN_FILAMENT_BINDING",                                                                       
"GOMF_MONOCARBOXYLIC_ACID_BINDING",                                                                  
"GOMF_UBIQUITIN_BINDING",                                                                            
"GOMF_ANTIOXIDANT_ACTIVITY",                                                                         
"GOMF_FATTY_ACID_BINDING",                                                                           
"GOMF_DNA_HELICASE_ACTIVITY",                                                                        
"GOMF_CYTOSKELETAL_PROTEIN_BINDING",                                                                 
"GOMF_SCAFFOLD_PROTEIN_BINDING",                                                                     
"GOMF_UBIQUITIN_LIKE_PROTEIN_TRANSFERASE_ACTIVITY",                                                  
"GOMF_ENDOPEPTIDASE_REGULATOR_ACTIVITY",                                                             
"GOMF_ENZYME_REGULATOR_ACTIVITY",                                                                    
"GOMF_PROTEIN_DIMERIZATION_ACTIVITY",                                                                
"GOMF_UBIQUITIN_LIKE_PROTEIN_BINDING",                                                               
"GOMF_PHOSPHOPROTEIN_PHOSPHATASE_ACTIVITY",                                                          
"GOMF_MAGNESIUM_ION_BINDING",                                                                        
"GOMF_ATP_DEPENDENT_ACTIVITY_ACTING_ON_DNA",                                                         
"GOMF_PHOSPHATASE_ACTIVITY",                                                                         
"GOMF_CLATHRIN_BINDING",                                                                             
"GOMF_ENZYME_INHIBITOR_ACTIVITY",                                                                    
"GOMF_PHOSPHATIDYLINOSITOL_BINDING",                                                                 
"GOMF_PROTEIN_TYROSINE_KINASE_BINDING",                                                              
"GOMF_TRANSCRIPTION_COREGULATOR_BINDING",                                                            
"GOMF_PHOSPHATIDYLINOSITOL_BISPHOSPHATE_BINDING",                                                    
"GOMF_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_CH_OH_GROUP_OF_DONORS",                                      
"GOMF_PHOSPHORIC_ESTER_HYDROLASE_ACTIVITY",                                                          
"GOMF_ENZYME_ACTIVATOR_ACTIVITY",                                                                    
"GOMF_CALCIUM_DEPENDENT_PROTEIN_BINDING",                                                            
"GOMF_HELICASE_ACTIVITY",                                                                            
"GOMF_PROTEIN_PHOSPHATASE_BINDING",                                                                  
"GOMF_DISORDERED_DOMAIN_SPECIFIC_BINDING", 
"GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME",                             
"GOMF_MRNA_3_UTR_BINDING",                                             
"GOMF_RRNA_BINDING",                                                   
"GOMF_OXIDOREDUCTION_DRIVEN_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
"GOMF_MRNA_BINDING",                                                   
"GOMF_PRE_MRNA_BINDING",                                               
"GOMF_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",                      
"GOMF_NUCLEOSOME_BINDING",                                             
"GOMF_PRIMARY_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",              
"GOMF_SINGLE_STRANDED_RNA_BINDING",                                    
"GOMF_PROTON_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",                      
"GOMF_ACTIVE_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY")) 

plot_NES(filtered_df = GOMF_combined, plottitle = "Ramos et al (2022) data - Control vs Schizophrenia Organoids - GO Molecular Functions")  



#################### January 17 2023

#Venn diagram of overlapping proteins between Ramos et al (2022) control/schizophrenia organoids and Sofia's control/FABP7 treated organoids

#read in the Ramos et al control and schizophrenia organoids dataset

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/schizophrenia/pathway analysis/")
ramos_ctr_scz_organoids_data = read.table("organoids.gct.txt", stringsAsFactors = FALSE, quote = "", header = TRUE, sep = "\t", check.names = FALSE)

setwd("C:/Users/qt09n/Desktop/Technical Analyst I UHN May 4 2021/organoid group/Sofia/pathway analysis/results")
sofia_ctr_fabp7<- read.table("FABP7_vehicle_control_data_processed.txt", stringsAsFactors = FALSE, sep = "\t", quote = "")

colnames(ramos_ctr_scz_organoids_data)<-ramos_ctr_scz_organoids_data[1,]

ramos_ctr_scz_organoids_data[,2]<- NULL

sofia_ctr_fabp7$gene <- rownames(sofia_ctr_fabp7)
colnames(ramos_ctr_scz_organoids_data)[1]<- "gene"
  
#find the number of common proteins between the Ramos dataset and Sofia's control/FABP7-treated organoids dataset
ctr_fabp7_proteins <- rownames(sofia_ctr_fabp7)
ctr_scz_organoid_proteins <- ramos_ctr_scz_organoids_data$NAME

common_proteins<- Reduce(intersect, list(ctr_scz_organoid_proteins, ctr_fabp7_proteins))

length(common_proteins)
#1635

#remove NA genes from the protein lists
ctr_scz_organoid_proteins<- na.omit(ctr_scz_organoid_proteins)
ctr_fabp7_proteins<-na.omit(ctr_fabp7_proteins)

#find the proteins which are unique to Sofia's control/FABP7 treated organoids dataset
sofia_ctr_fabp7_unique <- unique.comparisons(dataframe1 = sofia_ctr_fabp7, dataframe2 = ramos_ctr_scz_organoids_data)

#see which proteins are unique to Ramos et al (2022) control dataset compared to Sofia's control/FABP7 treated organoids
ramos_ctr_scz_organoids_data_unique <- unique.comparisons(dataframe1 = ramos_ctr_scz_organoids_data, dataframe2 = sofia_ctr_fabp7)

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
myCol <-c("#B3E2CD", "#FDCDAC")

# Generate 3 sets of 200 words
set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")

# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)



venn.diagram(
  x = list(ctr_fabp7_proteins, ctr_scz_organoid_proteins),
  category.names = c(paste0("Melliou Proteins Dataset - ", "\n", "Control & FABP7-treated", "\n", "Organoids"), paste0("Ramos et al (2022)", "\n", "Proteins Dataset", "\n", "Control &", "\n",  "Schizophrenia", "\n", "Organoids")),
  filename = "Proteins Control, FABP7-treated and Schizophrenia Organoids.png",
  output = TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-30, 30),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)



venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


### re-do Venn diagrams excluding control samples - inlcude only FABP7 and Schizophrenia organoids

sofia_ctr <- sofia_ctr_fabp7[,1:5]
ramos_ctr <- ramos_ctr_scz_organoids_data[,c(1, 15:22)]
