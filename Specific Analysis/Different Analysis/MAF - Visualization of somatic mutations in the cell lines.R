

### Install maftools package

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")
library(maftools)

library(utils)
library(tidyverse)
library(data.table)

library(maftools)


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

write.table(MutationsT, file = "MutationsT.csv", row.names = F, sep = "\t")

##
laml <- read.maf(maf ="C:/GitHub/project-02-group-05/MutationsT.csv", useAll = T, verbose = T)


laml


## Documents
getSampleSummary(laml)


getGeneSummary(laml)


write.mafSummary(maf = laml, basename = 'laml')

## Visualization


#Ploting MAF summary

plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(maf = laml, top = 15)
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





