################################# Gene Ontology analysis with FC biomarkers ##########################################

# load biomarkers over FC 
# matrix we will work with: 
biomarkers_FC_genes= row.names(biomarkers_FC)
biomarkers_FC_values = as.matrix(biomarkers_FC_values)

# load the library 
library(clusterProfiler)

#################################################################################################################

### translate the gene names 

library(org.Hs.eg.db)

# do the translation 
# index 2 because already done analysis with biomarkers from p.values 

translated.genes2= bitr(biomarkers_FC_genes,fromType="SYMBOL", toType="ENTREZID",OrgDb = org.Hs.eg.db) 
head(translated.genes2)

#################################################################################################################

#### Gene Ontology Classification
# Notes: BP for Biological Process, MF for Molecular Function, and CC for Cellular Component

library(DOSE)

# take the needed gene ID
gene2=translated.genes2$ENTREZID
head(gene2)

### do classification

## 1. Biological Process
ggo2.1 <- groupGO(gene= gene2, OrgDb = org.Hs.eg.db,  ont = "BP", level = 3, readable = FALSE)
head(summary(ggo2.1))
# visualization 
barplot(ggo2.1, drop=TRUE, showCategory=12)

## 2. Molecular Function
ggo2.2 <- groupGO(gene= gene2, OrgDb = org.Hs.eg.db,  ont = "MF", level = 3, readable = FALSE)
head(summary(ggo2.2))
# visualization 
barplot(ggo2.2, drop=TRUE, showCategory=12)

## 3. Cellular Component
ggo2.3 <- groupGO(gene= gene2, OrgDb = org.Hs.eg.db,  ont = "CC", level = 3, readable = FALSE)
head(summary(ggo2.3))
# visualization 
barplot(ggo2.3, drop=TRUE, showCategory=12)

#################################################################################################################

### enrich GO 

# only works with gene id ENSEMBL 

library(org.Hs.eg.db)

# create gene data frame ans translate it 

gene2.2 <- row.names(biomarkers_FC_values)
gene.df2 <- bitr(gene2.2, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

##############################################################################################################
### Cellular Component
ego2.1 <- enrichGO(gene         = gene.df2$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

head(summary(ego2.1))

# visualization 
dotplot(ego2.1, showCategory=3)
cnetplot(ego2.1)
emapplot(ego2.1)

##############################################################################################################

### Biological Process
ego2.2 <- enrichGO(gene         = gene.df2$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

head(summary(ego2.2))

# visualization 
dotplot(ego2.2, showCategory=10)
cnetplot(ego2.2)
emapplot(ego2.2)

##############################################################################################################

### Molecular Function
ego2.3 <- enrichGO(gene         = gene.df2$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

head(summary(ego2.3))

# visualization 
dotplot(ego2.3, showCategory=10)
cnetplot(ego2.3)
emapplot(ego2.3)
