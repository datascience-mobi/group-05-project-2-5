
library(readr)
library(rstudioapi)


wd = dirname(rstudioapi::getSourceEditorContext()$path)

#Load data:

  Untreated = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
  Treated = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
  Basal = readRDS(paste0(wd,"/data/CCLE_basalexpression.rds"))
  Copynumber = readRDS(paste0(wd,"/data/CCLE_copynumber.rds"))
  Mutations = readRDS(paste0(wd,"/data/CCLE_mutations.rds"))
  Sensitivity = readRDS(paste0(wd,"/data/NegLogGI50.rds"))
  Drug_annotation = read_tsv(paste0(wd,"/data/drug_annotation.tsv"))
  Cellline_annotation = read_tsv(paste0(wd,"/data/cellline_annotation.tsv"))
  Metadata = read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))
  ## /data/ to find the file in the directory!
  
########################## BROAD ANALYSIS ###############################################

##########################HIGHLY VARIABLE GENES########################################################

  
##Find genes which vary the most in Untreated and Treated celllines:
## Seurat packages need to be installed: install.packages('Seurat')
##For the seurat library to be available, digest packages need to be installed: install.packages("digest")
library(Seurat)
#BroadAnU = Broad analysis untreated

BroadAnU <- CreateSeuratObject(counts = Untreated, project = "BroadAnU")
BroadAnU <- FindVariableFeatures(BroadAnU, selection.method = "vst", nfeatures= 2000)

# Identify the 10 most highly variable genes
top10Untreated <- head(VariableFeatures(BroadAnU), 10)


    # plot variable features with and without labels
    plotUntreated <- VariableFeaturePlot(BroadAnU)
    plotUntreatedlabeled <- LabelPoints(plot = plot1, points = top10Untreated, repel = TRUE, xnudge = 0, ynudge = 0)


#The same with treated data:
BroadAnT <- CreateSeuratObject(counts = Treated, project = "BroadAnT")
BroadAnT <- FindVariableFeatures(BroadAnT, selection.method = "vst", nfeatures= 2000)


# Identify the 10 most highly variable genes:
top10Treated <- head(VariableFeatures(BroadAnT), 10)


    # plot variable features with and without labels
    plotTreated <- VariableFeaturePlot(BroadAnT)
    plotTreatedlabeled <- LabelPoints(plot = plot1, points = top10Treated, repel = TRUE, xnudge = 0, ynudge = 0)


## Can we find matches between Treated and Untreated regarding their most variable gene expression?
top10Untreated %in% top10Treated
# There are 7 matches! 

    # Which genes are in the Untreated but not in Treated
    setdiff(top10Untreated, top10Treated)  #  "CD24"  "KRT18" "SPARC"

    # Which genes are in the Treated but not in Untreated
    setdiff(top10Treated, top10Untreated)  #  "CRISP3"   "RGS1"   "LOC101928635///ALDH1A2"

# Can we say that for the matching genes, the treatment did not change the high variation?
###################################################################################################

