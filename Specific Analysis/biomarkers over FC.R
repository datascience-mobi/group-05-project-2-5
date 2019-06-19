############################## biomarkers over FC ###############################################################
#after compairison of biomarkers over p.Values and FC we decided not to trust them, because a low p.Value was 
# not correlated with a high FC
# This is because for the statistical t.test low values in data frame also give low p.Values which seem to be 
# significant 

# Because we want to use the genes with the biggest diffrence in gene expression after durg treatment, we will 
# now focus on finding them over the FC values 

#################################################################################################################
 
### load Creat Vorinostat data ! 

#################################################################################################################

### creat FC data 

FC <- TreatedVorinostat - UntreatedVorinostat

# work with mean of the rows because we only want to compare the genes 
FC_meanrow= rowMeans(FC)

#################################################################################################################

### sort the data 

# work with absolute value to find the highest values
# because we want to have the most up and down regulated genes 
FC_abs= abs(FC_meanrow)

## sort the values to get the 100 largest values 

sortedFC_abs <- sort(FC_abs, decreasing = TRUE)
sortedFC_abs <- as.matrix(sortedFC_abs)

# take the first 100 for biomarkers 
biomarkers_FC = sortedFC_abs[1:100,]
biomarkers_FC <- as.matrix(biomarkers_FC)

# see that the last ones have very similar values 

# creat vector with gene names 
biomarkers_FC_genes= row.names(biomarkers_FC)


#################################################################################################################

### creat matrix with FC values (positive and negativ) 

# add the abs values to FC matrix 
FC_both= cbind(FC_meanrow,FC_abs)
FC_both=as.data.frame(FC_both)

# order this matrix 
FC_both_sorted <- FC_both[order(FC_both$FC_abs, decreasing = TRUE),]

# FC values of biomarkers
# take the first 100 of the sorted matrix, should be the same as un biomarkers_FC 
biomarkers_FC_values = FC_both_sorted[1:100,]
# remove the absolute values
biomarkers_FC_values <- subset( biomarkers_FC_values, select = -FC_abs)
biomarkers_FC_values = as.matrix(biomarkers_FC_values)


#################################################################################################################

### Creat Matirx with biomarkers with FC of all samples 
FCbiomarkers <- FC[biomarkers_FC_genes,]

