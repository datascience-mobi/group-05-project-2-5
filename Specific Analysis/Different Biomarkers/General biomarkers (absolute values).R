####Find top100 general biomarkers for all celllines:

library(readr)
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)

#Load needed data:

Untreated_notnormalized = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
Treated_notnormalized = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
Metadata = read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))

#We need a FC matrix that contains only the data from Vorinostat-treated celllines:
UntreatedVorinostatcolumns <- grep(pattern = "vorinostat",colnames(Untreated_notnormalized))
TreatedVorinostatcolumns <- grep(pattern = "vorinostat",colnames(Treated_notnormalized))

UntreatedVorinostat <- Untreated_notnormalized[,UntreatedVorinostatcolumns]
TreatedVorinostat <- Treated_notnormalized[,TreatedVorinostatcolumns]

FCVorinostat <- TreatedVorinostat - UntreatedVorinostat


# work with mean of the rows because we only want to compare the genes:
FCVorinostatabs= abs(FCVorinostat)
FCVorinostatmean <- apply(FCVorinostatabs, 1, mean)

# sort the values to get the 100 largest values:
sortedgeneralbiomarker <- sort(FCVorinostatmean, decreasing = TRUE)
sortedgeneralbiomarker <- as.matrix(sortedgeneralbiomarker)

#select the top 100 general biomarkers:
top100generalbiomarkers = sortedgeneralbiomarker[1:100,]
top100generalbiomarkers <- as.matrix(top100generalbiomarkers)

#create vector with gene names:
generalbiomarkergenes = row.names(top100generalbiomarkers)

############################################################################################################

### compare this biomarkers to the biomarkers without absolute values 
 
### creat the other biomarkers 
FCVorinostatmean_noabs= rowMeans(FCVorinostat)

# work with absolute value to find the highest values
# because we want to have the most up and down regulated genes 
FC_abs= abs(FCVorinostatmean_noabs)

## sort the values to get the 100 largest values 

sortedFC_abs <- sort(FC_abs, decreasing = TRUE)
sortedFC_abs <- as.matrix(sortedFC_abs)

# take the first 100 for biomarkers 
biomarkers_FC = sortedFC_abs[1:100,]
biomarkers_FC <- as.matrix(biomarkers_FC)

# see that the last ones have very similar values 

# creat vector with gene names 
biomarkers_FC_genes= row.names(biomarkers_FC)

### compairison 

setequal(biomarkers_FC_genes,generalbiomarkergenes)
# returs false

diff1= setdiff(biomarkers_FC_genes,generalbiomarkergenes)
length(diff1)
# only 4 diffrent biomarkers 

diff1
# biomarkers in biomarkers_FC_genes but not in generalbiomarkergenes

diff2= setdiff(generalbiomarkergenes,biomarkers_FC_genes)
length(diff2)
# only 4 diffrent biomarkers 

diff2
# biomarkers in biomarkers_FC_genes but not in generalbiomarkergenes

