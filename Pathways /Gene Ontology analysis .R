
# work according to this paper: 
# https://bioconductor.statistik.tu-dortmund.de/packages/3.1/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.pdf

###################################################################################################
###################################################################################################

#install cluster profiler 
'if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")'

library(clusterProfiler)

###################################################################################################
###################################################################################################

############################## creat a list of genes/biomarkers ###################################

# pValues for biomarkers
pValues <- apply(VorinostatTotal, 1, function(x) t.test(x[col_untreated],x[col_treated],paired = TRUE, alternative = "two.sided")$p.value)

# sort the p.values 
sortedpValues <- sort(pValues, decreasing = FALSE)
sortedpValues <- as.matrix(sortedpValues)

# take the first 100 p.values for biomarkers 
biomarkers <- sortedpValues[1:100,]
biomarkers <- as.matrix(biomarkers)

# creat vector with genes 
biomarkers.genes = row.names(biomarkers)
biomarkers.genes

############################ translating amog diffrent gene ID types ##############################

# install needed libary form translation 
'if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")'

library(org.Hs.eg.db)

# do the translation 
translated.genes= bitr(biomarkers.genes,fromType="SYMBOL", toType="ENTREZID",OrgDb = org.Hs.eg.db) 
head(translated.genes)

# in which types the gene names can be translated by this libary 
idType("org.Hs.eg.db")

###################################################################################################
###################################################################################################

################################# Gene Ontology analysis ##########################################

###################################################################################################

################################### Gene Ontology Classification ##################################
# einheitliche Bezeichnungen in der Bioinformatik 

# install needed libary 
'if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DOSE")'

library(DOSE)
data(biomarkers)
gene=translated.genes$ENTREZID
head(gene)

# gene classification 
ggo <- groupGO(gene= gene, OrgDb = org.Hs.eg.db,  ont = "BP", level = 3, readable = FALSE)
ggo2=as.data.frame(ggo)               
head(summary(ggo2))
plot (ggo2)

# visualization 
barplot(ggo, drop=TRUE, showCategory=12)



'# transfrom GO name to symbol 
x=ggo2$ID
translated.ggo= bitr(x,fromType="GO", toType="SYMBOL",OrgDb = org.Hs.eg.db)
ggo.sym.genes=translated.ggo$SYMBOL
# bind to ggo matrix 
ggo.sym=cbind(ggo2,ggo.sym.genes)
# will return the symbol gene name of all GO genes, but we do not have a count for every GO gene'



# To Do: GO -> Pathway/Function, Symbol -> GO -> pathway  

##############################################################################################################
##############################################################################################################

########################################## enrich GO ########################################################
# Notes: BP for Biological Process, MF for Molecular Function, and CC for Cellular Component
# only works with gene id ENSEMBL, do not know why 

library(org.Hs.eg.db)
gene2 <- row.names(biomarkers)
gene.df <- bitr(gene2, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

##############################################################################################################
### Cellular Component
ego1 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

head(summary(ego1))

# visualization 
dotplot(ego1, showCategory=3)
cnetplot(ego1)
##############################################################################################################

### Biological Process
ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

head(summary(ego2))

# visualization 
dotplot(ego2, showCategory=10)
cnetplot(ego2)

##############################################################################################################

### Molecular Function
ego3 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

head(summary(ego3))

# visualization 
dotplot(ego3, showCategory=100)
cnetplot(ego3)









