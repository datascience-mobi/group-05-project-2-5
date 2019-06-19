# KEGG over-representation test 
###################################################################################################
# instrall needed package 
'if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGG.db")'

library(KEGG.db)

# see for other applications of this libary: 
# https://bioconductor.org/packages/release/data/annotation/manuals/KEGG.db/man/KEGG.db.pdf

###################################################################################################
# works better with biomarkers from FC
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

# problem: only finds 5 categories 