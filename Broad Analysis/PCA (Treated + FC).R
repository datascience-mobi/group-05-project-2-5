
library(readr)
library(rstudioapi)

#load data:

wd = dirname(rstudioapi::getSourceEditorContext()$path)

  Untreated = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
  Treated = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
  Basal = readRDS(paste0(wd,"/data/CCLE_basalexpression.rds"))
  Copynumber = readRDS(paste0(wd,"/data/CCLE_copynumber.rds"))
  Mutations = readRDS(paste0(wd,"/data/CCLE_mutations.rds"))
  Sensitivity = readRDS(paste0(wd,"/data/NegLogGI50.rds"))
  Drug_annotation = read_tsv(paste0(wd,"/data/drug_annotation.tsv"))
  Cellline_annotation = read_tsv(paste0(wd,"/data/cellline_annotation.tsv"))
  Metadata = read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))

  
########################## BROAD ANALYSIS ###############################################

##########################PRINCIPLE COMPONENT ANALYSIS###################################



#########PCA OVER TREATED DATA############################
  
treated.pca = prcomp(Treated, center=T, scale. = T)

  #use drug information from Metadata to color points in PCA:
  
        ##First we need to find out, which column in Metadata represents the different samples:
         samplecolumn <- (grep(colnames(Metadata), pattern = "sample"))
         samplecolumn
              #column 1 includes the samplenames
    
         
        #next, we need to check, if the sample-column includes the same celllines in the same order as the Treated matrix:
         dim(Metadata)
          #-> we have 2*819 rows since Metadata consists information for treated and untreated celllines.
    
      Metadatasamples <- Metadata[1:819,samplecolumn]
      identical(colnames(Treated), Metadatasamples)
       #--> TRUE, consequently the drug information of the Metadata-matrix can be assigned to the samples in the Treated-matrix sequentially:
  
  #Find the drug-column in Metadata:
  drugcolumn <- (grep(colnames(Metadata), pattern= "drug"))
  drugcolumn
  #column number 3 corresponds to drugs
      #for better readability, assign the column of interest to the name"Metadatadrugs"
       Metadatadrugs <- Metadata[1:819,drugcolumn]
 
      #Add Metadatadrugs as a new row to the Treated-matrix:
      Treatedwithdrugs <- rbind(Treated, Metadatadrugs)
        
         #save drug information as factors so it can be used for coloring:
         drugfactor <- as.factor(Treatedwithdrugs["Metadatadrugs",])
  
#####PLOT PCA & COLOR ACCORDING TO DRUG#####################################
        
     #15 different drugs, consequently we need 15 different colors:
      palette(rainbow(15))
        
        #Plot PCA: Principal component 1 and 2:
        plot(treated.pca$rotation[, 1], treated.pca$rotation[, 2], pch = 19, xlab = "PC1",ylab = "PC2", col=drugfactor, main = "Treated samples")
          
           #add a legend:
           levels <- as.factor(levels(drugfactor))
           legend("topright", inset = c(-0.4,0),levels(drugfactor), xpd = TRUE, pch=19, col = levels)
           par(mar=c(5, 4, 5, 8))
           
           
############PLOT PCA & COLOR ACCORDING TO TISSUE#########################
           
       #Find tissue-column in Metadata:
           tissuecolumn <- grep(colnames(Metadata), pattern= "tissue")
          Metadatatissue <- Metadata[1:819,tissuecolumn]
        
           #bind Metadatatissue as a new row to the Treated matrix:
           Treatedwithtissue <- rbind(Treated, Metadatatissue)
           
           #save tissue information as factors so it can be used for coloring:
           tissuefactor <- as.factor(Treatedwithtissue["Metadatatissue",])
           
           #since we have 9 different tissue types we need 9 different colors:
           palette(rainbow(9))
    
           
      ##plot PCA: PC 3 and PC 4 :
      plot(treated.pca$rotation[, 3], treated.pca$rotation[, 4], pch = 19, xlab = "PC1",ylab = "PC2", col=tissuefactor, main= "Treated samples")

            
           #add a legend:
                levels <- as.factor(levels(tissuefactor))
                legend("topright", inset = c(-0.3,0), levels(tissuefactor), xpd = TRUE, pch=19, col = levels) 
                 par(mar=c(5, 4, 5, 8))
           
                 
                 
                 
                 
###################################################################################################

          

#########PCA WITH THE FC MATRIX##########################################

FC <- Treated-Untreated

pca = prcomp(FC, center = T, scale. = T)
      #to get information about pca data: print(pca)
        #see how much variance is exolained by each principle component:
         plot(pca, type = "l")
      

  #Where are the vorinostat-treated celllines in the PCA plot? --> use ifelse-function
      Metadata <-as.data.frame(Metadata)    
      Marking <- ifelse(Metadata$drug == "vorinostat", "yellow", "black")
      #add the information, whether samples belong to Vorinostat, to the FC matrix:
      HighlightVorinostat <- cbind(`FC` = Marking)
          
        #plot the PCA:
          plot(pca$rotation[, 3], pca$rotation[, 4], col = HighlightVorinostat, pch = 19, xlab = "PC3", ylab = "PC4")

          
###########PLOT PCA & COLOR ACCORDING TO TISSUE###########################
  
      FCwithtissue <- rbind(FC, Metadatatissue)
      
          #save tissue information as factors so it can used for coloring:
           tissuefactor <- as.factor(FCwithtissue["Metadatatissue",])
         
           #9 different tissue types:
           palette(rainbow(9))
           
   #plot different PCs and see which combination groups the samples best:
    plot(pca$rotation[, 1], pca$rotation[, 2], col = tissuefactor, pch = 19, xlab = "PC1", ylab = "PC2")
    plot(pca$rotation[, 2], pca$rotation[, 3],col= tissuefactor, pch = 19, xlab = "PC2", ylab = "PC3")
    plot(pca$rotation[, 3], pca$rotation[, 4],col= tissuefactor, pch = 19, xlab = "PC3", ylab = "PC4")
    plot(pca$rotation[, 1], pca$rotation[, 3],col= tissuefactor, pch = 19, xlab = "PC1", ylab = "PC3")
    #-->different tissues do not seem to group the points in any PC combination:
 
    
#############PLOT PCA & COLOR ACCORDING TO DRUGS#########################
    
    FCwithdrugs <- rbind(FC, Metadatadrugs)
    
          #save drug information as factors so it can used for coloring:
          drugfactor <- as.factor(FCwithdrugs["Metadatadrugs",])

          #15 different tissue types:
          palette(rainbow(15))
          
    #adapt the size of the plot so that the legend can be seen:
    par(mar=c(5,4,5,9))
    
    #Plot PCA: PC 3 & 4, since they group different drugs best:
    plot(pca$rotation[, 3], pca$rotation[, 4], col = drugfactor , pch = 19, xlab = "PC3", ylab = "PC4", main = "FC PCA")

     
            #add a legend with the drugs:
            levels <- as.factor(levels(drugfactor))
             legend("topright", inset = c(-0.4,0), levels(drugfactor), xpd = TRUE, pch=19, col = levels)




