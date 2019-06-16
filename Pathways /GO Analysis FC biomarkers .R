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


