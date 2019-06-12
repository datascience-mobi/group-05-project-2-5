#T-test

#Create a matrix that  contains only the celllines treated with vorinostat:
# see Creat Vorinostat data ! 

##############################################################################################################
    #For the t-test, we need to check if genes and celllines are in the same order in both matrices:
    # kann auch mit 'oder' sortiert werden, falls nicht der Fall order function 
    identical(rownames(UntreatedVorinostat), rownames(TreatedVorinostat))
    #-> genes are in the same order
    
    identical(colnames(UntreatedVorinostat), colnames(TreatedVorinostat))
    #Colnames differ since the the matrices have a different drug dosis
    
      #Remove the information about the dosis (replace it with nothing):
      Untreatedcolnames <- sub(pattern = "_0nM_24h", "", colnames(UntreatedVorinostat))
      Treatedcolnames <- sub(pattern = "_5000nM_24h", "", colnames(TreatedVorinostat))
        #Check if the colnames are equal now:
       identical(Untreatedcolnames, Treatedcolnames)
       #-> since the colnames are equal now, matrices do not differ regarding their order of celllines, the dosis is the only difference, thus we can use both matrices to perform a t-test:
       
##############################################################################################################
              
 #T-tests require normal distributions. Since our sample is  very big, a normality test is not needed necessarily.
 #Anyway: Check normality of the FC-values via qqplot:
       FC <- TreatedVorinostat - UntreatedVorinostat
       qqnorm(FC)
       qqline(FC, col= "red")
       #We see, that we have a heavy tailed distribution. Nevertheless, we assume normality.
       

       #For normalization, a package needs to be installed: install.packages("BBmisc")
       library(BBmisc)
       FCnorm <- normalize(FC, method= "scale")
       #Check, if data look more normal distributed now:
       qqnorm(FCnorm)
       qqline(FCnorm, col= "red")
       #-> more linear than before

# same check for treated and untreated 
       qqnorm(TreatedVorinostat)
       qqline(TreatedVorinostat, col= "red")
       
       qqnorm(UntreatedVorinostat)
       qqline(UntreatedVorinostat, col= "red")

#Perform a paired t-test by using the apply-function:
       #H0 hypothesis: Gene expression does not change significantly after drug treatment.
       #H1 hypothesis: Gene expression changes significantly after drug treatment.
       
       dim(UntreatedVorinostat)
       #->59 celllines and 13299 genes
    
##############################################################################################################       
          
 #Create a common matrix for treated and untreated vorinostat data:
 VorinostatTotal <- cbind(UntreatedVorinostat, TreatedVorinostat)
  
  pValues1 <- apply(VorinostatTotal, 1, function(x) t.test(x[1:59], x[60:118],paired = TRUE, alternative = "two.sided")$p.value)

  
    # variable um Zahlen zu ersetzen 
 col_untreated = grep ('_0nM',colnames(VorinostatTotal))
 col_untreated
 
 col_treated = grep ('_5000nM',colnames(VorinostatTotal))
 col_treated

 
 ##############################################################################################################
 
 
  # gives a list with all information from t.test
  t.test.Vorinostat = apply(VorinostatTotal, 1, function(x) t.test(x[col_untreated], x[col_treated],paired = TRUE, alternative = "two.sided"))  
 
  pValues <- apply(VorinostatTotal, 1, function(x) t.test(x[col_untreated],x[col_treated],paired = TRUE, alternative = "two.sided")$p.value)
  sum(pValues < 0.05)   
  #-> gives a p value for each gene but takes the mean of rows
  
  #sort p-Values:
  sortedpValues <- sort(pValues, decreasing = FALSE)
  sortedpValues <- as.matrix(sortedpValues)
  SignificantsortedpValues[100,]
  
  #add p-Values as a new column to each gene:
  VorinostatwithpValues <- cbind(VorinostatTotal, pValues)
  dim(VorinostatwithpValues)
  
  
  ##############################################################################################################
  
  #select those rows with smallest p.Values:
  Biomarker <- VorinostatwithpValues[VorinostatwithpValues$pValues <=  sortedpValues[100,],]
  #-->Matrix with biomarker-genes!!
  
  
  #Sort biomarkers according to increasing p-Values (most significant p-Values on top of the data frame):
  Biomarker <- as.data.frame(Biomarker)
  Biomarkersorted <- Biomarker[order(Biomarker$pValues),]



  ##find 100 genes with smallest p-Values:
  #we want the smallest values at the beginning, thus they should be sorted by increasing value:
  sortedpValues <- sort(pValues, decreasing = FALSE)
  sortedpValues <- as.matrix(sortedpValues)
  
  #100 genes with the smallest p-Value are most significant and serve as our biomarkers:
  biomarkers <- sortedpValues[1:100,]
  #save biomarkers in a matrix;
  biomarkers <- as.matrix(biomarkers)
  
  
  
  '#-> gives a vector with the corresponding rownumbers in VorinostatTotal -matrix:
  i=1
  l <- vector("list", 100)
  while(i<=100) {
    l[[i]] <- (grep(rownames(VorinostatTotal), pattern= names(biomarkers[i,])))
    i=i+1
  }'
  
  ##############################################################################################################
  #Not important here, but for the sake of completeness:
  # t.test FC 
  t_v_FC= apply(FC,1,t.test)
  