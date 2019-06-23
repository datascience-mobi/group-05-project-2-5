
library(readr)
library(rstudioapi)

#load data:

wd = dirname(rstudioapi::getSourceEditorContext()$path)

  Untreated_notnormalized = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
  Treated_notnormalized = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
  Metadata = read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))

  
########################## BROAD ANALYSIS ###############################################

##########################PRINCIPLE COMPONENT ANALYSIS###################################

#--> normalize data:
  
  Untreated <- apply(Untreated_notnormalized, 2, function(x){
    (x - mean(x)) / sd(x)
  })
  
  
  Treated <- apply(Treated_notnormalized, 2, function(x){
    (x - mean(x)) / sd(x)
  })
  
  
  FC_notnormalized <- Treated_notnormalized - Untreated_notnormalized
  FC <- apply(FC_notnormalized, 2, function(x){
    (x - mean(x)) / sd(x)
  })
  

#########PCA OVER TREATED DATA############################
  
treated.pca = prcomp(Treated, center=T, scale. = T)

  #use drug information from Metadata to color points in PCA:
         
        #first we need to check, if the sample-column includes the same celllines in the same order as the Treated matrix:
         dim(Metadata)
          #-> we have 2*819 rows since Metadata consists information for treated and untreated celllines.
      
       #Check, if the sample order is equal in the Treated-matrix and in Metadata:
        Metadata <- as.data.frame(Metadata)
        Metadatasamples <- Metadata[1:819,"sample"]
        identical(colnames(Treated), Metadatasamples)
        #--> TRUE, consequently the drug information of the Metadata-matrix can be assigned to the samples in the Treated-matrix sequentially:
  
       Metadatadrugs <- Metadata[1:819,"drug"]
 
      #Add Metadatadrugs as a new row to the Treated-matrix:
      Treatedwithdrugs <- rbind(Treated, Metadatadrugs)
        
         #save drug information as factors so it can be used for coloring:
         drugfactor <- as.factor(Treatedwithdrugs["Metadatadrugs",])
  
         
#####PLOT PCA & COLOR ACCORDING TO DRUG#####################################
        
     #15 different drugs, consequently we need 15 different colors:
      palette(rainbow(15))
        
        #Plot PCA: Principal component 1 and 2:
        plot(treated.pca$rotation[, 1], treated.pca$rotation[, 2], pch = 19, xlab = "PC1",ylab = "PC2", col=drugfactor, main = "PCA Treated samples colored by drug")
          
           #add a legend:
           levels <- as.factor(levels(drugfactor))
           legend("topright", inset = c(-0.4,0),levels(drugfactor), xpd = TRUE, pch=19, col = levels, title = "Drugs")
           par(mar=c(5, 4, 5, 8))
           
           
############PLOT PCA & COLOR ACCORDING TO TISSUE#########################
           
       #Find tissue-column in Metadata:
 
          Metadatatissue <- Metadata[1:819,"tissue"]
        
           #bind Metadatatissue as a new row to the Treated matrix:
           Treatedwithtissue <- rbind(Treated, Metadatatissue)
           
           #save tissue information as factors so it can be used for coloring:
           tissuefactor <- as.factor(Treatedwithtissue["Metadatatissue",])
           
           #since we have 9 different tissue types we need 9 different colors:
           palette(rainbow(9))
    
           
      ##plot PCA: PC 2 and PC 3 :
      plot(treated.pca$rotation[, 2], treated.pca$rotation[, 3], pch = 19, xlab = "PC2",ylab = "PC3", col=tissuefactor, main= "PCA Treated samples colored by tissue")

            
           #add a legend:
                levels <- as.factor(levels(tissuefactor))
                legend("topright", inset = c(-0.3,0), levels(tissuefactor), xpd = TRUE, pch=19, col = levels, title = "Tissue") 
                 par(mar=c(5, 4, 5, 8))
           
                 
                 
                 
                 
###################################################################################################

          

#########PCA WITH THE FC MATRIX##########################################


pca.FC = prcomp(FC, center = T, scale. = T)
      #to get information about pca data: print(pca)
        #see how much variance is exolained by each principle component:
         plot(pca.FC, type = "l")
      

  #Where are the vorinostat-treated celllines in the PCA plot? --> use ifelse-function
      Metadata <-as.data.frame(Metadata)    
      Marking <- ifelse(Metadata$drug == "vorinostat", "yellow", "black")
         #add the information, whether samples belong to Vorinostat, to the FC matrix:
         HighlightVorinostat <- cbind(`FC` = Marking)
          
        #plot the PCA:
          plot(pca.FC$rotation[, 3], pca.FC$rotation[, 4], col = HighlightVorinostat, pch = 19, xlab = "PC3", ylab = "PC4", main = "Highlighted Vorinostat samples in FC PCA")

          #adda a legend: 
          legend("topright", inset = c(-0.3,0), legend = c("Vorinostat","Other drugs"), xpd = TRUE, pch=19, col = c("yellow", "black")) 
          par(mar=c(5, 4, 5, 8))
          
          
###########PLOT PCA & COLOR ACCORDING TO TISSUE###########################
  
      FCwithtissue <- rbind(FC, Metadatatissue)
      
          #save tissue information as factors so it can used for coloring:
           tissuefactor <- as.factor(FCwithtissue["Metadatatissue",])
         
           #9 different tissue types:
           palette(rainbow(9))
           
   #plot different PCs and see which combination groups the samples best:
    plot(pca.FC$rotation[, 1], pca.FC$rotation[, 2], col = tissuefactor, pch = 19, xlab = "PC1", ylab = "PC2")
    plot(pca.FC$rotation[, 2], pca.FC$rotation[, 3],col= tissuefactor, pch = 19, xlab = "PC2", ylab = "PC3")
    plot(pca.FC$rotation[, 3], pca.FC$rotation[, 4],col= tissuefactor, pch = 19, xlab = "PC3", ylab = "PC4")
    plot(pca.FC$rotation[, 1], pca.FC$rotation[, 3],col= tissuefactor, pch = 19, xlab = "PC1", ylab = "PC3")
    #-->different tissues do not seem to group the points in any PC combination:
 
    
#############PLOT PCA & COLOR ACCORDING TO DRUGS#########################
    
    FCwithdrugs <- rbind(FC, Metadatadrugs)
    
          #save drug information as factors so it can used for coloring:
          drugfactor <- as.factor(FCwithdrugs["Metadatadrugs",])

          #15 different tissue types:
          palette(rainbow(15))
          
    #adapt the size of the plot so that the legend can be seen:
    par(mar=c(5,4,5,9))
    
    #Plot different PCAs:
    plot(pca.FC$rotation[, 1], pca.FC$rotation[, 2], col = drugfactor , pch = 19, xlab = "PC1", ylab = "PC2", main = "PCA FC colored by drugs")
    plot(pca.FC$rotation[, 2], pca.FC$rotation[, 3], col = drugfactor , pch = 19, xlab = "PC2", ylab = "PC3", main = "PCA FC colored by drugs")
    plot(pca.FC$rotation[, 3], pca.FC$rotation[, 4], col = drugfactor , pch = 19, xlab = "PC3", ylab = "PC4", main = "PCA FC colored by drugs")

     
            #add a legend with the drugs:
            levels <- as.factor(levels(drugfactor))
             legend("topright", inset = c(-0.4,0), levels(drugfactor), xpd = TRUE, pch=19, col = levels, title = "Drugs")




