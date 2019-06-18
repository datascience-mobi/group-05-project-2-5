library(readr)
library(rstudioapi)

#load data:

wd = dirname(rstudioapi::getSourceEditorContext()$path)

Untreated = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
Treated = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
Drug_annotation = read_tsv(paste0(wd,"/data/drug_annotation.tsv"))
Cellline_annotation = read_tsv(paste0(wd,"/data/cellline_annotation.tsv"))
Metadata = read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))

##First, load Vorinoistat data from file: "Create Vorinostat data"
#This code goes on with FCnorm:

#Perform a one-sample t.test with normalized FC data:
T.testFCnorm <- apply(FCnorm, 1, t.test)

#Save the results in a data frame:
TResults <- lapply(T.testFCnorm, function(.tres) {      
  
  data.frame(t.value=.tres[1],dfs=.tres[2],conf.int1=.tres$conf.int[1],conf.int2=
               .tres$conf.int[2],p.value=.tres[3])
})

Tfinalresults <- do.call(rbind, TResults)

#For the vulcano plot, a package needs to be installed:

#if (!requireNamespace('BiocManager', quietly = TRUE))
  #install.packages('BiocManager')
#BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)


#Since we are interested in general biomarkers, we use the mean FC of the columns:
log2FCnorm <- apply(FCnorm, 1, mean)

#add the FCvalues as a new column to the data frame containing the results of the t.test:
FCnormandpvalues <- cbind(Tfinalresults, log2FCnorm)

#Create a Vulcanoplot:
EnhancedVolcano(FCnormandpvalues,
                lab = rownames(FCnormandpvalues),
                x = 'log2FCnorm',
                y = 'p.value',
                xlim = c(-5, 8))

#If we do not use the normalized FC values, the plot looks different:

log2FC <- apply(FC, 1, mean)

#add the FCvalues as a new column to the data frame containing the results of the t.test:
FCandpvalues <- cbind(Tfinalresults, log2FC)



keyvals <- rep('black', nrow(FCandpvalues))
names(keyvals) <- rep('FC<Top100', nrow(FCandpvalues))

keyvals[which(FCandpvalues$log2FC >= min(biomarkers_FC))] <- 'gold'
names(keyvals)[which(FCandpvalues$log2FC >= min(biomarkers_FC))] <- 'positive FC in Top100'

keyvals[which(FCandpvalues$log2FC <= -min(biomarkers_FC))] <- 'green'
names(keyvals)[which(FCandpvalues$log2FC <= -min(biomarkers_FC))] <- 'negative FC in Top100'

#Create a Vulcanoplot:
EnhancedVolcano(FCandpvalues,
                lab = rownames(FCandpvalues),
                x = 'log2FC',
                y = 'p.value',
                title = 'Significance versus fold change of all genes',
                selectLab = biomarkers_FC_genes,
                transcriptLabSize = 1.8,
                colOverride = keyvals,
                pCutoff = 10e-14,
                FCcutoff = min(biomarkers_FC)
              )



