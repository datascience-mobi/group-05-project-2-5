# KEGG over-representation test 
###################################################################################################
# instrall needed package 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGG.db")

library(KEGG.db)

# see for other applications of this libary: 
# https://bioconductor.org/packages/release/data/annotation/manuals/KEGG.db/man/KEGG.db.pdf

###################################################################################################

kk <- enrichKEGG(gene         = translated.genes,
                 organism     = "human",
                 pvalueCutoff = 0.05,
                 readable= TRUE, 
                 use_internal_data = TRUE)

# here also problem with gene Id 