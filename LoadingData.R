####imports von libraries##############################################
library(readr)
#######################################################################


####Pfad setzen: jeder dort, wo die Daten auf dem PC liegen########################################################
setwd('~/Desktop/Uni/4.Semester/projekt/data/')
print('Aktuelles Verzeichniss: ') getwd()
#getwd() um auszugeben wo R momentan ausgef??hrt wird, bzw von wo der relative path gestartet wird
#######################################################################


##################Einlesen der ganzen Datens??tze#######################

#### Drug perturbation data from cell lines (transcriptome modulation and inhibition of cell growth)
##Hauptpaper
# Gene expression
mat_NCI_TPW_gep_untreated <- readRDS("NCI_TPW_gep_untreated.rds")
mat_NCI_TPW_gep_treated <- readRDS("NCI_TPW_gep_treated.rds")
tsv_NCI_TPW_metadata <- read_tsv("NCI_TPW_metadata.tsv")
# Drug sensitivity assay
mat_NegLogGI50.rds <- readRDS("NegLogGI50.rds")

####Basal molecular profiles of cancer cell lines 
##(mutations + gene copy number alterations + basal gene expression):
##nature paper
mat_CCLE_basalexpression <- readRDS("CCLE_basalexpression.rds")
mat_CCLE_copynumber.rds <- readRDS("CCLE_copynumber.rds")
mat_CCLE_mutations.rds <- readRDS("CCLE_mutations.rds")

#### Feature annotation
# cell line metadata
tsv_cellline_annotation <- read_tsv("cellline_annotation.tsv")
# mechanism of action, ect.
tsv_drug_annotation <- read_tsv("drug_annotation.tsv")


#######################################################################

#daten umwandeln
mat_NCI_TPW_gep_treated <- as.data.frame(mat_NCI_TPW_gep_treated)
# um Datentyp zu testen: is.recursive(mat_NCI_TPW_gep_treated)

mat_NCI_TPW_gep_untreated <- as.data.frame(mat_NCI_TPW_gep_untreated)

mat_NegLogGI50.rds <- as.data.frame(mat_NegLogGI50.rds)


