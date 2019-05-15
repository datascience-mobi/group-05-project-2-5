####imports von libraries##############################################
library(readr)
library(rstudioapi)
#######################################################################


#### find local directory #####
wd = dirname(rstudioapi::getSourceEditorContext()$path)
#########################


##################Einlesen der ganzen Daten#######################

#### Drug perturbation data from cell lines (transcriptome modulation and inhibition of cell growth)
##Hauptpaper
# Gene expression
 
NCI_TPW_gep_untreated <- readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
NCI_TPW_gep_treated <- readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
NCI_TPW_metadata <- read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))
# Drug sensitivity assay
NegLogGI50 <- readRDS(paste0(wd,"/data/NegLogGI50.rds"))

####Basal molecular profiles of cancer cell lines 
##(mutations + gene copy number alterations + basal gene expression):
##nature paper
CCLE_basalexpression <- readRDS(paste0(wd,"/data/CCLE_basalexpression.rds"))
CCLE_copynumber <- readRDS(paste0(wd,"/data/CCLE_copynumber.rds"))
CCLE_mutations <- readRDS(paste0(wd,"/data/CCLE_mutations.rds"))

#### Feature annotation
# cell line metadata
cellline_annotation <- read_tsv(paste0(wd,"/data/cellline_annotation.tsv"))
# mechanism of action, ect.
drug_annotation <- read_tsv(paste0(wd,"/data/drug_annotation.tsv"))


#######################################################################

#daten umwandeln
NCI_TPW_gep_treated <- as.data.frame(NCI_TPW_gep_treated)
# um Datentyp zu testen: is.recursive(mat_NCI_TPW_gep_treated)
NCI_TPW_gep_untreated <- as.data.frame(NCI_TPW_gep_untreated)

NegLogGI50<- as.data.frame(NegLogGI50)





