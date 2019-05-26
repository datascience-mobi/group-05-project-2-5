
library(readr)
library(rstudioapi)


wd = dirname(rstudioapi::getSourceEditorContext()$path)

  Untreated = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
  Treated = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
  Basal = readRDS(paste0(wd,"/data/CCLE_basalexpression.rds"))
  Copynumber = readRDS(paste0(wd,"/data/CCLE_copynumber.rds"))
  Mutations = readRDS(paste0(wd,"/data/CCLE_mutations.rds"))
  Sensitivity = readRDS(paste0(wd,"/data/NegLogGI50.rds"))
  Drug_annotation = read_tsv(paste0(wd,"/data/drug_annotation.tsv"))
  Cellline_annotation = read_tsv(paste0(wd,"/data/cellline_annotation.tsv"))
  Metadata = read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))
  ## /data/ to find the file in the directory!
  
########################## BROAD ANALYSIS ###############################################

##################################################################################################

########## colored with "level function" --> see loading data
############### Do we have batches ? ########
boxplot(Treated,medcol="red", border = NA, col= Metadata$drug, 
        xlab="drugs", ylab="gene expression", main= "gene expression over treated cellines", las=3)
  
#### Ja, da wir nicht ueberall einen gleichen Median haben und eindeutig Boxen erkennen, die jeweils 
# zu einem Medikament geh??ren 
  
  
########## remove batches ##############
#### load limma package  befor: http://bioconductor.org/packages/release/bioc/html/limma.html
library(limma)
#### which parameters are needed??? -> does not work yet 
removeBatchEffect()

##################################################################################################

##################################################################################################

##Find genes which vary the most in Untreated and Treated celllines:
## Seurat packages need to be installed: install.packages('Seurat')
library(Seurat)

#BroadAnU = Broad analysis untreated

BroadAnU <- CreateSeuratObject(counts = Untreated, project = "BroadAnU")
BroadAnU <- FindVariableFeatures(BroadAnU, selection.method = "vst", nfeatures= 2000)

# Identify the 10 most highly variable genes
top10Untreated <- head(VariableFeatures(BroadAnU), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(BroadAnU)
plot2 <- LabelPoints(plot = plot1, points = top10Untreated, repel = TRUE, xnudge = 0, ynudge = 0)
# Change from top1o to top10Untreated

#The same with treated data:
BroadAnT <- CreateSeuratObject(counts = Treated, project = "BroadAnT")
BroadAnT <- FindVariableFeatures(BroadAnT, selection.method = "vst", nfeatures= 2000)

# Identify the 10 most highly variable genes
top10Treated <- head(VariableFeatures(BroadAnT), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(BroadAnT)
plot2 <- LabelPoints(plot = plot1, points = top10Treated, repel = TRUE, xnudge = 0, ynudge = 0)
# Change from top10 to top10Treated

## Can we find matches between Treated and Untreated regarding their most variable gene expression?
top10Untreated == top10Treated

#There are matches! 4 genes are positioned equally in the top10 positions. A closer look at the top10 genes shows,that there are even 6 out of 10 genes which appear in both samples.

###################################################################################################

#####################################################################################################

###PCA with untreated celllines:

all.genes <- rownames(BroadAnU)
BroadAnU <- ScaleData(BroadAnU, features = all.genes)


BroadAnU <- RunPCA(BroadAnU, features = VariableFeatures(object = BroadAnU))
print(BroadAnU[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(BroadAnU, dims = 1:2, reduction = "pca")
DimPlot(BroadAnU, reduction = "pca")
DimHeatmap(BroadAnU, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(BroadAnU, dims = 2, cells = 500, balanced = TRUE)
ElbowPlot(BroadAnU)

#Problem: There is no clear elbow in the curve.

####################################################################################################

###################################################################################################

##### PCA over "normal" data ######
treated.pca = prcomp(Treated, center=T, scale. = T)
# colored drug
plot(treated.pca$rotation[, 1], treated.pca$rotation[, 2], pch = 19,
     xlab = "PC1",ylab = "PC2", col=Metadata$drug)
# colored tissue
plot(treated.pca$rotation[, 1], treated.pca$rotation[, 2], pch = 19,
     xlab = "PC1",ylab = "PC2", col=Metadata$tissue)
### --> found groups 

###################################################################################################


###################################################################################################

#Density plot which compares the mean expression of all celllines with and without treatment:
plot(density(colMeans(Untreated)) ,xlab = "Mean", main = "Effects of treatment on gene expression")
lines(density(colMeans(Treated)), col = "red")
legend(x = 1, y = 1, legend = c("Untreated", "Treated"), bty = "n", lty = 1, col = c("black", "red"), x.intersp = 0.5)



##Find those celllines which differ the most regarding their overall gene expression:

UntreatedT <- t(Untreated)
# calculate the distances and put the calculations into an object called distances
distances <- dist(UntreatedT)
# convert this distances object into a matrix. 
distances.m <- data.matrix(distances)

# we can extract the size of the object and the titles
dim <- ncol(distances.m)
names <- row.names(distances.m)
# now to create the visualisation of the difference matrix. 
# first the coloured boxes
image(1:dim, 1:dim, distances.m, axes = FALSE, xlab = "", ylab = "")
# now label the axis
axis(3, 1:dim, names, cex.axis = 0.8, las=3)
axis(2, 1:dim, names, cex.axis = 0.8, las=1)
# we could add the values of the differences with text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", distances.m), cex=1) but that is too confusing.


####################################################################################################################

####################################################################################################

### PCA with the foldchange matrix:

FC <- Treated-Untreated

pca = prcomp(FC, center = T, scale. = T)
      print(pca)
      plot(pca, type = "l")
      
          plot(pca$rotation[, 1], pca$rotation[, 2], pch = 19, xlab = "PC1", ylab = "PC2")
          plot(pca$rotation[, 2], pca$rotation[, 3], pch = 19, xlab = "PC2", ylab = "PC3")
          plot(pca$rotation[, 3], pca$rotation[, 4], pch = 19, xlab = "PC3", ylab = "PC4")
          #--> PC 3+4 seem to separate the points best.

#Where are the vorinostat-treated celllines in the PCA plot? --> use ifelse-function
cb1 <- ifelse(Metadata$drug == "vorinostat", "yellow", "black")
cb <- cbind(`FC` = cb1)
plot(pca$rotation[, 3], pca$rotation[, 4], col = cb, pch = 19, xlab = "PC3", ylab = "PC4")


#color the PCA points according to the cellline tissue:
tissue <- Metadata[1:819,6]
tissue <- t(tissue)
FC <- rbind(FC, tissue)
tissue <- as.factor(FC[13300,])

    plot(pca$rotation[, 1], pca$rotation[, 2], col = tissue, pch = 19, xlab = "PC1", ylab = "PC2")
    #--> I tried different PCs but could not find groups of cellines from the same tissue.

    # here also the levels could be used for coloring 
    # plot(pca$rotation[, 1], pca$rotation[, 2], col = Metadata$tissue, pch = 19, xlab = "PC1", ylab = "PC2")

##Separate points according to the 15 different drugs:
drug <- Metadata[1:819,3]
drug <- t(drug)
FC <- rbind(FC, drug)
drug <- as.factor(FC[13301,])

plot(pca$rotation[, 3], pca$rotation[, 4], col = drug, pch = 19, xlab = "PC3", ylab = "PC4")

  levels <- as.factor(levels(drug))
  legend("topright", inset = c(-0.3,0), levels(drug), xpd = TRUE, pch=19, col = levels )

#--> to see the legend I had to adapt the size of the plot with par(mar=c(3,4,3,7)+1)



