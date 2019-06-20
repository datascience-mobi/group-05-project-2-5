##PCA only using FC of biomarkers and Vorinostatsamples:


#Select only those genes from the FC matrix which were identified as biomarkers:
FCbiomarkers <- FC[biomarkers_FC_genes,]

#Execute the PCA:
FCbiomarkers.pca = prcomp(FCbiomarkers, center=T, scale. = T)

#For coloring according to tissue type we use information from Metadata:
Metadata <- as.data.frame(Metadata)

#We only need the tissueinformation of those samples treated with Vorinostat. Since Metadata includes them twice, once with and once without drug treatment, we select only the Treated samplenames:
Vorinostatsamples <- grep(Metadata$sample, pattern= "vorinostat_5000nM")
Metadatavorinostattissue <- Metadata[Vorinostatsamples,"tissue"]

#Add the tissueinformation as a new row to the Biomarkermatrix:
FCbiomarkerswithtissue <- rbind(FCbiomarkers, Metadatavorinostattissue)

#save tissue information as factors so it can be used for coloring:
vorinostattissuefactor <- as.factor(FCbiomarkerswithtissue["Metadatavorinostattissue",])

#9 different tissue types:
palette(rainbow(9))

#PC 3 & 4 seem to group the tissues best:
plot(FCbiomarkers.pca$rotation[, 3], FCbiomarkers.pca$rotation[, 4], col = vorinostattissuefactor, pch = 19, xlab = "PC3", ylab = "PC4", main = "PCA with FC of biomarkers colored according to tissue")

#Add a legend to the plot:
levels <- as.factor(levels(vorinostattissuefactor))
legend("topright", inset = c(-0.4,0),levels(Metadatavorinostattissue), xpd = TRUE, pch=19, col = levels)
par(mar=c(5, 4, 5, 9))
