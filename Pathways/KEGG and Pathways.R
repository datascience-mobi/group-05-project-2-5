# KEGG over-representation test 
###################################################################################################
# instrall needed package 
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("KEGG.db")

# see for other applications of this libary: 
# https://bioconductor.org/packages/release/data/annotation/manuals/KEGG.db/man/KEGG.db.pdf

###################################################################################################
library(KEGG.db)
library(clusterProfiler)
library(org.Hs.eg.db)
gene2 <- row.names(biomarkers_FC)
gene.df <- bitr(gene2, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene=gene.df$ENTREZID,pvalueCutoff = 0.05)
head(summary(kk))


# visualization
barplot(kk,showCategory=12)
dotplot(kk, showCategory=12)
cnetplot(kk,categorySize="geneNum")
emapplot(kk)


kk2 <- enrichMKEGG(gene=gene.df$ENTREZID,pvalueCutoff = 0.05)
head(summary(kk2))

# visualization
barplot(kk2,showCategory=12)
# only one categorie found 

###################################################################################################
# following enrichments based on this web page: 
# https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html#pathway-enrichment-analysis
###################################################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")

library(ReactomePA)

x <- enrichPathway(gene=gene.df$ENTREZID, pvalueCutoff = 0.05, readable = T )
head(as.data.frame(x))

# visualization
barplot(x,showCategory=12)
dotplot(x, showCategory=12)
cnetplot(x,categorySize="geneNum")
emapplot(x)


