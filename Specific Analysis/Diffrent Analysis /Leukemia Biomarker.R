####Find top100 biomarkers of Leukemia-celllines

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


#We need Metadata to see which of the celllines belong to Leukemia:

Metadata <- as.data.frame(Metadata)
#We only need the tissueinformation of those samples treated with Vorinostat. Since Metadata includes them twice, once with and once without drug treatment, we select only the Treated samplenames:
Vorinostatsamples <- grep(Metadata$sample, pattern= "vorinostat_5000nM")
Metadatavorinostattissue <- Metadata[Vorinostatsamples,"tissue"]

#Add the tissueinformation as a new row to the FCVorinostat-matrix:
FCVorinostatwithtissue <- rbind(FCVorinostat, Metadatavorinostattissue)

#select only those samples/celllines which belong to Leukemia:
FCVorinostatwithtissue <- as.data.frame(FCVorinostatwithtissue)
Leukemiasamples <- colnames(FCVorinostatwithtissue[Metadatavorinostattissue== "Leukemia"])
FCVorinostatLeukemia <- FCVorinostatwithtissue[,Leukemiasamples]

#Now we do not need the tissue information any longer:
FCVorinostatLeukemia <- FCVorinostatLeukemia[-13300,]

#We want the FC-mean for each gene:
FCVorinostatLeukemiaabs <- abs(FCVorinostatLeukemia)
FCVorinostatLeukemiameans <- apply(FCVorinostatLeukemiaabs,1 ,mean)

#Sort the biomarkers that those with the biggest FC are on top:
sortedFCleukemiamean <- sort(FCVorinostatLeukemiameans, decreasing = TRUE)
sortedFCleukemiamean <- as.matrix(sortedFCleukemiamean)

# take the first 100 as Leukemia-biomarkers:
biomarkersLeukemia = sortedFCleukemiamean[1:100,]
biomarkersLeukemia <- as.matrix(biomarkersLeukemia)

# create a vector with gene names of Leukemia-biomarkers: 
biomarkersLeukemiaGenes = row.names(biomarkersLeukemia)



#Comparison with general up- and downregulating biomarkers:
intersect(biomarkersLeukemiaGenes, biomarkers_FC_genes)

#Comparison with general absolute value biomarkers:
intersect(biomarkersLeukemiaGenes, generalbiomarkergenes)
biomarkersLeukemiaGenes[biomarkersLeukemiaGenes == generalbiomarkergenes]

