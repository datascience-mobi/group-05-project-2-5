########################### Vorinostat target genes in FC matrix ###############################################

################################################################################################################
library(readr)
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)

Drug_annotation = read.table(paste0(wd,"/data/drug_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)

################################################################################################################

# search for HDAC in biomarkers 
target <- grep(pattern = "HDAC",colnames(top100generalbiomarkers))
length(target)

# In drug annotation HDAC genes are mentioned as target of vorinostat, but these genes are not in our defined 
# biomarkers, so know we want to check their FC values 

###find target of vorinostat 

target=Drug_annotation$target
target=as.data.frame(target)

# take row with vorinostat 
target_vorinostat=target[9,]
target_vorinostat

# new vector with names as strings, because otherwise problem with spacer | 
target_vorinostat=c("HDAC1","HDAC10","HDAC11","HDAC2","HDAC3","HDAC5","HDAC6","HDAC8","HDAC9")


### find FC values of targets 

# creat FC data 
FC =TreatedVorinostat - UntreatedVorinostat
FC= rowMeans(FC)
FC=as.data.frame(FC)
genes=row.names(FC)
FC_new=cbind(FC,genes)

# filter data 
library(dplyr)
FC_target=as.data.frame(filter(FC_new, genes %in% target_vorinostat))

FC_target
# the change in gene expression is really low, most of the genes are up regulated insted of being down regulated
# maybe more production of HDAC because they are blocked by vorinostat? 




