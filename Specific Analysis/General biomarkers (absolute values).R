####Find top100 general biomarkers for all celllines:

library(readr)
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)

#Load needed data:

Untreated = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
Treated = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
Metadata = read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))

#We need a FC matrix that contains only the data from Vorinostat-treated celllines:
UntreatedVorinostatcolumns <- grep(pattern = "vorinostat",colnames(Untreated))
TreatedVorinostatcolumns <- grep(pattern = "vorinostat",colnames(Treated))

UntreatedVorinostat <- Untreated[,UntreatedVorinostatcolumns]
TreatedVorinostat <- Treated[,TreatedVorinostatcolumns]

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

