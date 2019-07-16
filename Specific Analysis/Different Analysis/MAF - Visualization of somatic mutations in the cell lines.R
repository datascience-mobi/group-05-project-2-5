

### Install maftools package

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")

library(maftools)


### Install other necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("VariantAnnotation")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

## Load libraries

library(utils)
library(tidyverse)
library(data.table)

library(Biostrings)
library(maftools)
library(ComplexHeatmap)

## Set working directory

wd = ("C:/GitHub/project-02-group-05/")
getwd()=='C:/GitHub/project-02-group-05'


## Edit the file so that it can be read by the package

Mutations <- readRDS(paste0(wd,"/data/CCLE_mutations.rds"))

names(Mutations)[names(Mutations) == "Tumor_Seq_Allele1"] <- "Tumor_Seq_Allele2"
names(Mutations)[names(Mutations) == "Start_position"] <- "Start_Position"
names(Mutations)[names(Mutations) == "End_position"] <- "End_Position"
names(Mutations)[names(Mutations) == "Hugo_Symbol"] <- "Hugo_Symbol"
rownames(Mutations) <- c()
MutationsT <- Mutations[,c(1,4,5,6,10,11,8,9,16,2,3,7,12,13,14,15,17,18)]


##Create a file with all genes
write.table(MutationsT, file = "MutationsT.csv", row.names = F, sep = "\t")

##
laml <- read.maf(maf ="C:/GitHub/project-02-group-05/MutationsT.csv", useAll = T, verbose = T)


laml


## Create a File with 100 Biomarkers
BM_mut = MutationsT[ which((MutationsT$Hugo_Symbol) 
                                  %in% rownames(biomarkers_FC_values100)), ]

write.table(BM_mut, file = "BM_mut.csv", row.names = F, sep = "\t")

BM_laml <- read.maf(maf ="C:/GitHub/project-02-group-05/BM_mut.csv", useAll = T, verbose = T)


#####ALL MUTATIONS ####


## Documents
getSampleSummary(laml)


getGeneSummary(laml)


write.mafSummary(maf = laml, basename = 'laml')

## Visualization


#Ploting MAF summary

plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

## Ploting oncoplot with Transversions/Transitions

oncoplot(maf = laml, top = 15, draw_titv = TRUE)

#oncoplot(maf = laml, top = 15, draw_titv = TRUE, additionalFeature = c("Tumor_Sample_Barcode", "MALME-3M"))

#Altered in all samples as this are all cancer cell lines

##Oncoplot expression values

set.seed(seed = 1024)
exprs_tbl = data.frame(genes = getGeneSummary(x = laml)[1:20, Hugo_Symbol],
                       exprn = rnorm(n = 10, mean = 12, sd = 5))
head(exprs_tbl)

oncoplot(maf = laml, exprsTbl = exprs_tbl)

#CHANGE EXPRESSION VALUES


## Oncostrips for the top10 biomarkers
oncostrip(maf = laml, genes = c('DHRS2',
                                'ABAT', 
                                'SERPINI1', 
                                'MIR612///NEAT1	',
                                'HBA2///HBA1	',
                                'CLU',
                                'NMI',
                                'STC1',
                                'AREG',
                                'NSMAF',
                                'SERPINH1')
          )



## Oncostrips for the 10 genes with the least change in expression



oncostrip(maf = laml, genes = c('ADAMTS6', 
                                'YWHAQ',
                                'EMID1', 
                                'PLCH2', 
                                'HRH2',
                                'WNT10B',
                                'RGPD6///RGPD8///RGPD3///RGPD4///RGPD5',
                                'ANGPTL8',
                                'CABP2',
                                'UPK2',
                                'CLPS')
        )

## Transitions and Transversion

laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)

##Mutation load vs TCGA cohorts
laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML')


## Somatic interactions
somaticInteractions(maf = laml, top = 30, pvalue = c(0.05, 0.1))

summary(somaticInteractions(maf = laml, top = 30, pvalue = c(0.05, 0.1)))

pairlist(somaticInteractions(maf = laml, top = 30, pvalue = c(0.05, 0.1)))

##                gene_set     pvalue
##1: LRP1B, PCDH15, MT-ND5 0.04964611

oncostrip(maf = laml, genes = c('LRP1B', 'PCDH15', 'MT-ND5'))









###### ONLY BM ####

BM_laml

## Documents
getSampleSummary(BM_laml)


getGeneSummary(BM_laml)


write.mafSummary(maf = laml, basename = 'laml')

## Visualization


#Ploting MAF summary

plotmafSummary(maf = BM_laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

## Ploting oncoplot

oncoplot(maf = BM_laml, top = 15)
#Altered in all samples as this are all cancer cell lines


## Oncostrips for the top10 biomarkers
oncostrip(maf = laml, genes = c('DHRS2',
                                'ABAT', 
                                'SERPINI1', 
                                'MIR612///NEAT1	',
                                'HBA2///HBA1	',
                                'CLU',
                                'NMI',
                                'STC1',
                                'AREG',
                                'NSMAF',
                                'SERPINH1')
          )


## Transitions and Transversion

BM.laml.titv = titv(maf = BM_laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = BM.laml.titv)

##Mutation load vs TCGA cohorts
laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML')


# Drug Gene interactions

dgi = drugInteractions(maf = laml, fontSize = 0.75)

dnmt3a.dgi = drugInteractions(genes = "MUC16", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]


OncogenicPathways(maf = laml)

PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS")

PlotOncogenicPathways(maf = laml, pathways = "NOTCH")

PlotOncogenicPathways(maf = laml, pathways = "TGF-Beta")


### Mutational signatures

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")


library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocInstaller::biocValid()

biocLite("stringi")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")



library(BiocManager)

library(BSgenome.Hsapiens.UCSC.hg19)



laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
plotApobecDiff(tnm = laml.tnm, maf = laml, pVal = 0.2)



