library(readr)
library(rstudioapi)


wd = dirname(rstudioapi::getSourceEditorContext()$path)

  #Load data

    Untreated = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
    Treated = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
    Basal = readRDS(paste0(wd,"/data/CCLE_basalexpression.rds"))
    Copynumber = readRDS(paste0(wd,"/data/CCLE_copynumber.rds"))
    Mutations = readRDS(paste0(wd,"/data/CCLE_mutations.rds"))
    Sensitivity = readRDS(paste0(wd,"/data/NegLogGI50.rds"))

    Metadata = read.table(paste0(wd,"/data/NCI_TPW_metadata.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
    Cellline_annotation = read.table(paste0(wd,"/data/cellline_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
    Drug_annotation = read.table(paste0(wd,"/data/drug_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
    
###BOXPLOT 
  #Treated celllines:
    
    par(mar= c(5, 4, 4, 7) + 0.1)
    drug <- Metadata$drug
    palette(rainbow(15))
    boxplot(Treated, medcol="black", border = drug, col= drug, xlab="drugs", ylab="gene expression", main= "Gene expression treated celllines", names= FALSE, xaxt= "n", boxwex=1, boxlty =0)
    
      #add a legend to see which color corresponds to which drug:
       levels <- as.factor(levels(drug))
       legend("topright", inset = c(-0.3,0), legend= levels(drug), xpd = TRUE, pch=19, col = levels, title = "drugs")
  
  #Boxplot for untreated celllines:
    par(mar= c(5, 4, 4, 7) + 0.1)
    drug <- Metadata$drug
    palette(rainbow(15))
    boxplot(Untreated, medcol="black", border = drug, col= drug, xlab="drugs", ylab="gene expression", main= "Gene expression untreated celllines", names= FALSE, xaxt= "n", boxwex=1, boxlty =0)
    
      levels <- as.factor(levels(drug))
      legend("topright", inset = c(-0.3,0), legend= levels(drug), xpd = TRUE, pch=19, col = levels, title = "drugs")
  
  #-> big differences in gene expression between different drug groups before drug treatment confirm that we have batches.

    