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
 
Untreated <- readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
Treated <- readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
#NCI_TPW_metadata <- read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))
Metadata = read.table(paste0(wd,"/data/NCI_TPW_metadata.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
# Drug sensitivity assay
Sensitivity <- readRDS(paste0(wd,"/data/NegLogGI50.rds"))

####Basal molecular profiles of cancer cell lines 
##(mutations + gene copy number alterations + basal gene expression):
##nature paper
Basal <- readRDS(paste0(wd,"/data/CCLE_basalexpression.rds"))
Copynumber <- readRDS(paste0(wd,"/data/CCLE_copynumber.rds"))
Mutations <- readRDS(paste0(wd,"/data/CCLE_mutations.rds"))

#### Feature annotation
# cell line metadata
#cellline_annotation <- read_tsv(paste0(wd,"/data/cellline_annotation.tsv"))
Cellline_annotation = read.table(paste0(wd,"/data/cellline_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
# mechanism of action, ect.
#drug_annotation <- read_tsv(paste0(wd,"/data/drug_annotation.tsv"))
Drug_annotation = read.table(paste0(wd,"/data/drug_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)


#######################################################################

#daten umwandeln
Treated <- as.data.frame(Treated)
# um Datentyp zu testen: is.recursive(mat_NCI_TPW_gep_treated)
Untreated <- as.data.frame(Untreated)

Sensitivity<- as.data.frame(Sensitivity)

################### auf levels zugreifen: levels(datensatz$spalte)
levels(Metadata$drug)





