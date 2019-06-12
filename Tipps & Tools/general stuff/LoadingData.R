####imports von libraries##############################################
library(readr)
library(rstudioapi)
#######################################################################


#### find local directory #####
wd = dirname(rstudioapi::getSourceEditorContext()$path)
#########################


##################Einlesen der ganzen Daten#######################
 
Untreated <- readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
Treated <- readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))

Metadata = read.table(paste0(wd,"/data/NCI_TPW_metadata.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)

Sensitivity <- readRDS(paste0(wd,"/data/NegLogGI50.rds"))


Basal <- readRDS(paste0(wd,"/data/CCLE_basalexpression.rds"))
Copynumber <- readRDS(paste0(wd,"/data/CCLE_copynumber.rds"))
Mutations <- readRDS(paste0(wd,"/data/CCLE_mutations.rds"))

Cellline_annotation = read.table(paste0(wd,"/data/cellline_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
Drug_annotation = read.table(paste0(wd,"/data/drug_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)


#######################################################################

#daten umwandeln
Treated <- as.data.frame(Treated)
# um Datentyp zu testen: is.recursive(mat_NCI_TPW_gep_treated)
Untreated <- as.data.frame(Untreated)

Sensitivity<- as.data.frame(Sensitivity)

################### auf levels zugreifen: levels(datensatz$spalte)
levels(Metadata$drug)






