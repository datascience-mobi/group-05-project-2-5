install.packages("GSEABase")
library(GSEABase)

BiocManager::install(c("genefilter"))

library(Biobase)
library(genefilter)

library(limma)
library(RColorBrewer)
library(GSVA)

cacheDir <- system.file("extdata", package="GSVA") 
cachePrefix <- "cache4vignette_"
file.remove(paste(cacheDir, list.files(cacheDir, pattern=cachePrefix), sep="/"))




gsva_document<- gsva(data.matrix(FC),
                      generalbiomarkergenes,
                      method="gsva",
                      mx.diff=TRUE,
                      verbose=TRUE,
                      parallel.sz=8)

#Set colour and breaks for heatmap
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
myBreaks <- seq(-1.5, 1.5, length.out=101)

heat <- t(scale(t(gsva_document)))

heatmap(heat,
        col=myCol,
        breaks=myBreaks,
        main="Title",
        key=TRUE,
        keysize=1.0,
        key.title="",
        key.xlab="Enrichment Z-score",
        scale="none", 
        density.info="none",
        reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
        trace="none",
        cexRow=1.0,
        cexCol=1.0,
        distfun=function(x) dist(x, method="euclidean"),
        hclustfun=function(x) hclust(x, method="ward.D2"))