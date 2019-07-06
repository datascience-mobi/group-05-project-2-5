if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVAdata")

generalbiomarkergenes2=as.list(generalbiomarkergenes2)
FC=TreatedVorinostat-UntreatedVorinostat

library(GSVA)
g<- gsva(FC,
         generalbiomarkergenes2,
         mx.diff=TRUE,
         verbose=TRUE)



library(GSVAdata)
data(c2BroadSets)
g<- gsva(FC,
         generalbiomarkergenes2,
     c2BroadSets,
     mx.diff=TRUE,
     verbose=TRUE)

library(GSVAdata)
data(leukemia)
g<- gsva(FC,
         generalbiomarkergenes2,
         leukemia,
         mx.diff=TRUE,
         verbose=TRUE)




heat <- t(scale(t(g)))

heatmap(heat,
        col=myCol,
        breaks=myBreaks,
        main="GSVA Enrichment Score of our Biomarker",
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


library(gplots)
heatmap.2(heat,
          col=myCol,
          breaks=myBreaks,
          main="GSVA Enrichment Score of our Biomarker",
          labRow="", labCol="",
          xlab="genes",
          ylab="samples",
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



