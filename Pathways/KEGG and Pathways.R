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
top100generalbiomarkers=as.data.frame(top100generalbiomarkers)
colnames(top100generalbiomarkers)[1] <- "fold Change"
top100generalbiomarkers=cbind(top100generalbiomarkers,generalbiomarkergenes)
colnames(top100generalbiomarkers)[2] <- "genes"
gene2 <- top100generalbiomarkers$genes
gene.df <- bitr(gene2, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene=gene.df$ENTREZID, keyType="kegg", pvalueCutoff = 0.05)
head(summary(kk))


# visualization
barplot(kk,showCategory=12)
dotplot(kk, showCategory=12)
cnetplot(kk,categorySize="geneNum")
emapplot(kk)


kk2 <- enrichMKEGG(gene=gene.df$ENTREZID,pvalueCutoff = 0.05, pAdjustMethod = "none")
head(summary(kk2))
# visualization
barplot(kk2,showCategory=12)
# only one categorie found 

# finde more categories 
kk3 <- enrichKEGG(gene=gene.df$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "none")
head(summary(kk3))
barplot(kk3,showCategory=100)
dotplot(kk3, showCategory=12)
cnetplot(kk3,categorySize="geneNum")
emapplot(kk3)



# colored FC 
library(enrichplot)
library(DOSE)
top100generalbiomarkers=as.data.frame(top100generalbiomarkers)
colnames(top100generalbiomarkers)[1] <- "fold Change"
top100generalbiomarkers=cbind(top100generalbiomarkers,generalbiomarkergenes)
colnames(top100generalbiomarkers)[2] <- "genes"

FC <- top100generalbiomarkers$`fold Change`
names(FC) <- gene.df$ENTREZID

cnetplot(kk3,categorySize="pvalue",foldChange=FC)
heatplot(kk3, foldChange=FC)
upsetplot(kk3)


## view pathway
'if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")'

library(pathview)
pathview(gene.data = FC, 
         pathway.id = "hsa05202", 
         species = "hsa", 
         limit = list(gene=5, cpd=1))



###################################################################################################
# following enrichments based on this web page: 
# https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html#pathway-enrichment-analysis
###################################################################################################
'if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")'

library(ReactomePA)

x <- enrichPathway(gene=gene.df$ENTREZID, pvalueCutoff = 0.05, readable = T )
head(as.data.frame(x))

# visualization
barplot(x,showCategory=12)
dotplot(x, showCategory=12)
cnetplot(x,categorySize="geneNum")
emapplot(x)




