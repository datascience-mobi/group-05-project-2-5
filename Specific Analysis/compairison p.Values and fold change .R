################################ compairison p.Values and fold change ###########################

#################################################################################################

#################################################################################################

### create data for p.values ###

pValues <- apply(VorinostatTotal, 1, function(x) t.test(x[col_untreated],x[col_treated],paired = TRUE, alternative = "two.sided")$p.value)

# sort the p.values 
sortedpValues <- sort(pValues, decreasing = FALSE)
sortedpValues <- as.matrix(sortedpValues)

# take the first 100 p.values for biomarkers 
biomarkers <- sortedpValues[1:100,]
biomarkers <- as.matrix(biomarkers)

# creat vector with genes 
biomarkers.genes = row.names(biomarkers)
# check if it has the correct length
length(biomarkers.genes)

#################################################################################################

#################################################################################################

### create data for FC ###

FC <- TreatedVorinostat - UntreatedVorinostat

# work with mean of the rows because we only want to compare the genes 
FC_meanrow= rowMeans(FC)

# sort the FC 
sortedFC <- sort(FC_meanrow, decreasing = FALSE)
sortedFC <- as.matrix(sortedFC)

# take the first 100 values like done above for compairison 
FC100<- sortedFC[1:100,]
FC100 <- as.matrix(FC100)

# creat vector with genes 
FC100.genes = row.names(FC100)
# check if it has the correct length
length(FC100.genes)

#################################################################################################

#################################################################################################

### compare genes vectors ### 

setequal(FC100.genes,biomarkers.genes)
# returs FALSE --> not the same genes 

diff= setdiff(FC100.genes,biomarkers.genes)
length(diff)
# returs 73 diffrent genes 

# diffrent biomarkers over p.Values from t.test and over FC values 

#################################################################################################

#################################################################################################

### compare p.Values and FC values in a new matrix ###

#creat new matrix
pV_FC= cbind(FC_meanrow,sortedpValues)

head(pV_FC)
# a high p.value does not lead to a high FC value

#################################################################################################

#################################################################################################

### creat diffrence matrix to see how big the differences are ### 

diff_pV_FC= as.matrix(FC_meanrow-sortedpValues)
summary(diff_pV_FC)

# Which genes differ most? 

# creat new data, because sort function will not return names 
rm(diff_pV_FC)

FC <- TreatedVorinostat - UntreatedVorinostat
FC_meanrow= rowMeans(FC)

pValues <- apply(VorinostatTotal, 1, function(x) t.test(x[col_untreated],x[col_treated],paired = TRUE, alternative = "two.sided")$p.value)
sortedpValues <- sort(pValues, decreasing = FALSE)

diff_pV_FC= (FC_meanrow-sortedpValues)

# sort the data, positve vales 
mostdiff_max= sort(diff_pV_FC, decreasing = T)[1:10]
mostdiff_max=as.matrix(mostdiff_max)

# sort the data, negativ vales 
mostdiff_min= sort(diff_pV_FC, decreasing = F)[1:10]
mostdiff_min=as.matrix(mostdiff_min)


# To Do: Wich differences between FC and p.Values do we accept?

