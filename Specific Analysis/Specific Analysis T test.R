#T-test

#Create a matrix that  contains only the celllines treated with vorinostat:

    
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
       
 #T-tests require normal distributions. Since our sample is  very big, a normality test is not needed necessarily.
 #Anyway: Check normality of the FC-values via qqplot:
       FC <- TreatedVorinostat - UntreatedVorinostat
       qqnorm(FC)
       qqline(FC, col= "red")
       #We see, that we have a heavy tailed distribution. Nevertheless, we assume normality.
       
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
       
 
 VorinostatTotal <- cbind(UntreatedVorinostat, TreatedVorinostat)
  
  # variable um Zahlen zu ersetzen 
 col_untreated = grep ('_0nM',colnames(VorinostatTotal))
 col_untreated
 
 col_treated = grep ('_5000nM',colnames(VorinostatTotal))
 col_treated

  
  pValues <- apply(VorinostatTotal, 1, function(x) t.test(x[col_untreated],x[col_treated],paired = TRUE, alternative = "two.sided")$p.value)
  sum(pValues < 0.05)   
  #Problem: gives a p value for each gene but takes the mean of rows
  
  t.test.Vorinostat = apply(VorinostatTotal, 1, function(x) t.test(x[1:59], x[60:118],paired = TRUE, alternative = "two.sided")) 
  # gives a list with all informations from t.test
  
  #MARGIN= c(1,2) means that the t-test should be performed over rows and columns
  i <- 1
  if (i<=59) {
    pValues(i) <- apply(VorinostatTotal, MARGIN = c(1,2), function(x) t.test(x= VorinostatTotal[,i], y= VorinostatTotal[,59+i],paired = TRUE, alternative = "two.sided")$p.value)
    i= i+1
  }
    #-> does not work, computer calculates for hours
  