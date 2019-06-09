
# work according to this paper: 
# https://bioconductor.statistik.tu-dortmund.de/packages/3.1/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.pdf

###################################################################################################
###################################################################################################

#install cluster profiler 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

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
###################################################################################################

################################### Gene Ontology Classification ##################################
# einheitliche Bezeichnungen in der Bioinformatik 

# install needed libary 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DOSE")

library(DOSE)

# gene classification 
ggo <- groupGO(gene= biomarkers.genes, OrgDb = org.Hs.eg.db,  ont = "BP", level = 3, readable = FALSE)
ggo2=as.data.frame(ggo)               
head(summary(ggo))
plot (ggo2)

# visualization 
barplot(ggo, drop=TRUE, showCategory=1)

## here we can see that something is wrong with the counts 


# transfrom GO name to symbol 
x=ggo2$ID
translated.ggo= bitr(x,fromType="GO", toType="SYMBOL",OrgDb = org.Hs.eg.db)
ggo.sym.genes=translated.ggo$SYMBOL

# bind to ggo matrix 
ggo.sym=cbind(ggo2,ggo.sym.genes)

# problem: not the same size, some GO have more then one gene 
# because GO are groups???

# To Do: GO -> Pathway/Function, Symbol -> GO -> pathway  

#################################### GO over-representation test ###################################

## hypergeometric model to assess wether a numer of selected genes is associated with disease is 
# lager than expected

# only works with ENTREZID Gene IdType 
# creat vector with ENTREZID Genes

# contains symbol and ENTREZID (see above)
translated.genes
# vector with only ENTREZID Genes
e.genes= translated.genes$ENTREZID
e.genes=as.data.frame(e.genes)
e.genes

# does not work ! :(
# problem with gene id 
# perform GO Enrichment Analysis 
ego <- enrichGO(e.genes,  OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC",
                pvalueCutoff = 0.01, pAdjustMethod = "BH", universe,
                qvalueCutoff = 0.05,readable = TRUE)

ego2=as.data.frame(ego)
head(summary(ego2))
                










