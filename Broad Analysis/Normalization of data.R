library(readr)
library(rstudioapi)


wd = dirname(rstudioapi::getSourceEditorContext()$path)

#Load data

Metadata = read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))
Untreated_notnormalized = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
Treated_notnormalized = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))


#Normalize data:
  
Untreated <- apply(Untreated_notnormalized, 2, function(x){
  (x - mean(x)) / sd(x)
})


Treated <- apply(Treated_notnormalized, 2, function(x){
  (x - mean(x)) / sd(x)
})


FC_notnormalized <- Treated_notnormalized - Untreated_notnormalized
FC <- apply(FC_notnormalized, 2, function(x){
  (x - mean(x)) / sd(x)
})



#Before normalization, we could not see drug groups when we plotted PC1 against PC2. Check, if FC PCA looks better now:

pca = prcomp(FC, center = T, scale. = T)

#####PLOT PCA & COLOR ACCORDING TO DRUG:

Metadata <- as.data.frame(Metadata)
Metadatadrugs <- Metadata[1:819,"drug"]
FCwithdrugs <- rbind(FC, Metadatadrugs)

drugfactor <- as.factor(FCwithdrugs["Metadatadrugs",])

palette(rainbow(15))
par(mar=c(5,4,5,9))

##Try PC 1 & 2:
plot(pca$rotation[, 1], pca$rotation[, 2], col = drugfactor , pch = 19, xlab = "PC1", ylab = "PC2", main = "FC PCA colored by drugs")
levels <- as.factor(levels(drugfactor))
legend("topright", inset = c(-0.4,0), levels(drugfactor), xpd = TRUE, pch=19, col = levels)


#Try PC 2 & 3:
plot(pca$rotation[, 2], pca$rotation[, 3], col = drugfactor , pch = 19, xlab = "PC1", ylab = "PC2", main = "FC PCA colored by drugs")
levels <- as.factor(levels(drugfactor))
legend("topright", inset = c(-0.4,0), levels(drugfactor), xpd = TRUE, pch=19, col = levels)

#---> after normalization we can see a more clear grouping of celllines which were treated with the same drug!
#We should use the normalized data for all further analyzes.
