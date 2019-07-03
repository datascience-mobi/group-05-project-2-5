if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVAdata")

library(GSVAdata)


library(clusterProfiler)
library(org.Hs.eg.db)
top100generalbiomarkers=as.data.frame(top100generalbiomarkers)
colnames(top100generalbiomarkers)[1] <- "fold Change"
top100generalbiomarkers=cbind(top100generalbiomarkers,generalbiomarkergenes)
colnames(top100generalbiomarkers)[2] <- "genes"
gene2 <- top100generalbiomarkers$genes
gene.df <- bitr(gene2, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
gene=gene.df$ENTREZID
gene=as.list(gene)

top100generalbiomarkers=cbind(top100generalbiomarkers,gene)
top100generalbiomarkers=as.matrix.ExpressionSet(top100generalbiomarkers)

FC= TreatedVorinostat-UntreatedVorinostat
generalbiomarkergenes2=as.list(generalbiomarkergenes)

data(c2BroadSets)
g<- gsva(top100generalbiomarkers,gene,
     c2BroadSets,
     mx.diff=TRUE,
     verbose=TRUE)



jjjf=gsva(FC,generalbiomarkergenes2,
          c2BroadSets,
     method=c("gsva", "ssgsea", "zscore", "plage"),
     kcdf=c("Gaussian", "Poisson", "none"),
     abs.ranking=FALSE,
     min.sz=1,
     max.sz=Inf,
     parallel.sz=0,
     parallel.type="SOCK",
     mx.diff=TRUE,
     tau=switch(method, gsva=1, ssgsea=0.25, NA),
     ssgsea.norm=TRUE,
     verbose=TRUE)


heat <- t(scale(t(g)))

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


library(gplots)
heatmap.2(heat,
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
