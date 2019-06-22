if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("enrichplot")

library(enrichplot)

library(KEGG.db)
library(clusterProfiler)
library(org.Hs.eg.db)
gene2 <- row.names(biomarkers_FC)
gene.df <- bitr(gene2, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene=gene.df$ENTREZID, keyType="kegg", pvalueCutoff = 0.05)
head(summary(kk))

cnetplot(kk,circular=TRUE)
heatplot(kk)
upsetplot(kk)

# finde more categories 
kk3 <- enrichKEGG(gene=gene.df$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "none")
cnetplot(kk3,categorySize="geneNum",circular=TRUE)
heatplot(kk3)
upsetplot(kk3)