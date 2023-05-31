#may 24 2023
# Ramos GSM6720583 17 pcw germinal matrix sample - differential expression/pathway analysis

#try to plot out the more interesting pathways from the Enrichr pathways table for 
#Gene Ontology Molecular Pathways (GOMF), Gene Ontology Biological Processes(GOBP),
#and Gene Ontology Cellular components (GOCC)
 

setwd("C:/Users/Diamandis Lab II/Documents/Queenie")

source("protein_organoid_working_matrix_main_analysis_functions.R")

library(enrichR)

#read in the results from Enrichr package 
setwd("C:/Users/Diamandis Lab II/Documents/Queenie/ramos/output/March 23 2023/GSM6720853/differential expression")

enriched.highFABP7 <-readRDS("enriched.highFABP7.rds")
enriched.lowFABP7 <- readRDS("enriched.lowFABP7.rds")

#GO Biological Processes pathways for high FABP7 Glutamatergic neurons
enriched.highFABP7.GOBP <- enriched.highFABP7[["GO_Biological_Process_2015"]]

enriched.lowFABP7.GOBP <- enriched.lowFABP7[["GO_Biological_Process_2015"]]

#do a comparisons test to see what pathways are unique to the high FABP7 glutamatergic neurons group
high_FABP7_GOBP_unique <-unique.comparisons.enrichr(enriched.highFABP7.GOBP, enriched.lowFABP7.GOBP)

#do a comparisons test to see what pathways are unique to the high FABP7 glutamatergic neurons group
low_FABP7_GOBP_unique <-unique.comparisons.enrichr(enriched.lowFABP7.GOBP, enriched.highFABP7.GOBP)

#find out pathways which are common to both high and low FABP7 glutamatergic neuron groups
common_pathways_high_and_low_fabp7_glutamatergic_neurons<-common.comparisons.enrichr(enriched.highFABP7.GOBP, enriched.lowFABP7.GOBP)


#pick out some of the more interesting pathways that show up in the results for plotting
highFABP7.GOBP.terms <- c("neuron projection guidance (GO:0097485)", "axon guidance (GO:0007411)",
                          "regulation of neuron projection development (GO:0010975)",
                          "regulation of actin polymerization or depolymerization (GO:0008064)",
                          "regulation of neuron differentiation (GO:0045664)",
                          "negative regulation of ERBB signaling pathway (GO:1901185)",
                          "positive regulation of protein complex assembly (GO:0031334)",
                          "cerebral cortex radially oriented cell migration (GO:0021799)",
                          "neurotrophin TRK receptor signaling pathway (GO:0048011)",
                          "actin cytoskeleton organization (GO:0030036)",
                          "dendrite morphogenesis (GO:0048813)",
                          "glial cell development (GO:0021782)",
                          "axon cargo transport (GO:0008088)",
                          "cerebral cortex cell migration (GO:0021795)",
                          "somatic cell DNA recombination (GO:0016444)",
                          "regulation of neuron apoptotic process (GO:0043523)",
                          "glial cell proliferation (GO:0014009)",
                          "glucocorticoid receptor signaling pathway (GO:0042921)",
                          "regulation of cell morphogenesis involved in differentiation (GO:0010769)",
                          "neuron differentiation (GO:0030182)",
                          "positive regulation of neurogenesis (GO:0050769)",
                          "neurofilament cytoskeleton organization (GO:0060052)",
                          "forebrain cell migration (GO:0021885)",
                          "anterograde axon cargo transport (GO:0008089)",
                          "negative regulation of neuron apoptotic process (GO:0043524)",
                          "regulation of oligodendrocyte differentiation (GO:0048713)",
                          "neural tube formation (GO:0001841)",
                          "cell morphogenesis involved in neuron differentiation (GO:0048667)",
                          "astrocyte development (GO:0014002)",
                          "positive regulation of neuron projection development (GO:0010976)",
                          "regulation of neuron death (GO:1901214)",
                          "neuron development (GO:0048666)",
                          "positive regulation of nervous system development (GO:0051962)",
                          "neuroepithelial cell differentiation (GO:0060563)",
                          "neuron recognition (GO:0008038)",
                          "positive regulation of neuron differentiation (GO:0045666)",
                          "fatty acid biosynthetic process (GO:0006633)",
                          "cell differentiation in spinal cord (GO:0021515)",
                          "central nervous system neuron differentiation (GO:0021953)",
                          "positive regulation of calcium ion transmembrane transporter activity (GO:1901021)",
                          "regulation of dendrite development (GO:0050773)",
                          "neutral lipid catabolic process (GO:0046461)",
                          "dendritic cell migration (GO:0036336)",
                          "regulation of gliogenesis (GO:0014013)",
                          "regulation of synaptic plasticity (GO:0048167)",
                          "spinal cord motor neuron differentiation (GO:0021522)",
                          "regulation of neuronal synaptic plasticity (GO:0048168)",
                          "regulation of glial cell differentiation (GO:0045685)",
                          "positive regulation of lipid kinase activity (GO:0090218)",
                          "protein localization to plasma membrane (GO:0072659)",
                          "positive regulation of transmembrane transport (GO:0034764)",
                          "regulation of axon extension (GO:0030516)",
                          "cellular response to fatty acid (GO:0071398)",
                          "cell projection morphogenesis (GO:0048858)",
                          "glycerolipid catabolic process (GO:0046503)",
                          "negative regulation of calcium ion-dependent exocytosis (GO:0045955)",
                          "synaptic vesicle maturation (GO:0016188)",
                          "regulation of low-density lipoprotein particle receptor biosynthetic process (GO:0045714)",
                          "substrate-independent telencephalic tangential interneuron migration (GO:0021843)",
                          "neurotransmitter secretion (GO:0007269)",
                          "neuron-neuron synaptic transmission (GO:0007270)",
                          "neuron migration (GO:0001764)",
                          "germinal center formation (GO:0002467)",
                          "positive regulation of astrocyte differentiation (GO:0048711)",
                          "neuron fate determination (GO:0048664)",
                          "positive regulation of potassium ion transmembrane transporter activity (GO:1901018)",
                          "regulation of neurotransmitter levels (GO:0001505)",
                          "regulation of dendritic spine development (GO:0060998)",
                          "regulation of calcineurin-NFAT signaling cascade (GO:0070884)",
                          "positive regulation of cAMP biosynthetic process (GO:0030819)",
                          "establishment of vesicle localization (GO:0051650)",
                          "positive regulation of fat cell differentiation (GO:0045600)",
                          "cellular response to calcium ion (GO:0071277)",
                          "positive regulation of oligodendrocyte differentiation (GO:0048714)",
                          "astrocyte differentiation (GO:0048708)",
                          "myelin maintenance (GO:0043217)",
                          "detection of calcium ion (GO:0005513)",
                          "neurotransmitter transport (GO:0006836)",
                          "neurogenesis (GO:0022008)",
                          "peripheral nervous system axon ensheathment (GO:0032292)",
                          "regulation of cholesterol storage (GO:0010885)",
                          "regulation of lipoprotein particle clearance (GO:0010984)",
                          "dendritic spine organization (GO:0097061)",
                          "myelination in peripheral nervous system (GO:0022011)",
                          "neurotransmitter uptake (GO:0001504)",
                          "regulation of synapse assembly (GO:0051963)",
                          "regulation of lipid kinase activity (GO:0043550)",
                          "cell proliferation in forebrain (GO:0021846)",
                          "protein localization to synapse (GO:0035418)",
                          "dopaminergic neuron differentiation (GO:0071542)",
                          "regulation of voltage-gated calcium channel activity (GO:1901385)",
                          "cerebral cortex development (GO:0021987)",
                          "neutral lipid biosynthetic process (GO:0046460)",
                          "positive regulation of sodium ion transmembrane transport (GO:1902307)",
                          "regulation of glial cell proliferation (GO:0060251)",
                          "forebrain neuron differentiation (GO:0021879)",
                          "central nervous system development (GO:0007417)",
                          "axonogenesis (GO:0007409)",
                          "axon regeneration (GO:0031103)",
                          "dendritic cell chemotaxis (GO:0002407)",
                          "cellular response to lipid (GO:0071396)",
                          "calcium-mediated signaling (GO:0019722)",
                          "neural precursor cell proliferation (GO:0061351)",
                          "regulation of synaptic transmission, glutamatergic (GO:0051966)",
                          "positive regulation of neuron apoptotic process (GO:0043525)",
                          "long-chain fatty-acyl-CoA biosynthetic process (GO:0035338)",
                          "transmission of nerve impulse (GO:0019226)",
                          "regulation of calcium ion transmembrane transport (GO:1903169)",
                          "axonal fasciculation (GO:0007413)",
                          "regulation of synaptic vesicle exocytosis (GO:2000300)",
                          "positive regulation of lipid storage (GO:0010884)",
                          "response to fatty acid (GO:0070542)",
                          "regulation of neurotransmitter secretion (GO:0046928)",
                          "positive regulation of axonogenesis (GO:0050772)",
                          "acetyl-CoA metabolic process (GO:0006084)",
                          "neutral lipid catabolic process (GO:0046461)",
                          "glycerolipid catabolic process (GO:0046503)",
                          "axon guidance (GO:0007411)",
                          "neuron projection guidance (GO:0097485)",
                          "regulation of cell projection organization (GO:0031344)"
                          )

enriched.highFABP7.GOBP.selected <- enriched.highFABP7.GOBP[enriched.highFABP7.GOBP$Term %in% highFABP7.GOBP.terms,]

#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html

plotEnrich(enriched.highFABP7.GOBP.selected, showTerms = 50, numChar = 50, y = "Count", orderBy = "P.value")

enriched.lowFABP7.GOBP

#pick out some of the more interesting pathways that show up in the results for plotting
#for the low FABP7 Glutamatergic Neurons Cluster group
enriched.lowFABP7.GOBP.terms <- (c("positive regulation of GTPase activity (GO:0043547)",
                                   "regulation of cell projection organization (GO:0031344)",
                                   "neuron projection guidance (GO:0097485)",
                                   "axon guidance (GO:0007411)",
                                   "synaptic transmission (GO:0007268)",
                                   "regulation of cell morphogenesis involved in differentiation (GO:0010769)",
                                   "cell projection morphogenesis (GO:0048858)",
                                   "regulation of ion transmembrane transport (GO:0034765)",
                                   "neuron recognition (GO:0008038)",
                                   "neuron projection morphogenesis (GO:0048812)",
                                   "calcium ion transmembrane transport (GO:0070588)",
                                   "regulation of transmembrane transport (GO:0034762)",
                                   "regulation of membrane potential (GO:0042391)",
                                   "regulation of axonogenesis (GO:0050770)",
                                   "axonogenesis (GO:0007409)",
                                   "positive regulation of neuron differentiation (GO:0045666)",
                                   "positive regulation of nervous system development (GO:0051962)",
                                   "positive regulation of neurogenesis (GO:0050769)",
                                   "calcium ion transport (GO:0006816)",
                                   "positive regulation of cell projection organization (GO:0031346)",
                                   "neurotrophin signaling pathway (GO:0038179)",
                                   "positive regulation of neuron projection development (GO:0010976)",
                                   "positive regulation of nervous system development (GO:0051962)",
                                   "regulation of synaptic transmission (GO:0050804)",
                                   "axonogenesis (GO:0007409)",
                                   "positive regulation of neuron projection development (GO:0010976)",
                                   "regulation of membrane potential (GO:0042391)",
                                   "positive regulation of neurogenesis (GO:0050769)",
                                   "regulation of cell-substrate adhesion (GO:0010810)",
                                   "regulation of cell junction assembly (GO:1901888)",
                                   "synapse organization (GO:0050808)"
                                   ))

enriched.lowFABP7.GOBP.selected <- enriched.lowFABP7.GOBP[enriched.lowFABP7.GOBP$Term %in% enriched.lowFABP7.GOBP.terms,]

plotEnrich(enriched.lowFABP7.GOBP.selected, showTerms = 60, numChar = 50, y = "Count", orderBy = "P.value")

