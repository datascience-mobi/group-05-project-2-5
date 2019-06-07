######################################## variance ########################################

### fold change data

# creat variance data
var_v_FC=(apply(FC,1,var))
var_v_FC

# sort the data
# important: do not store as.data.frame or as.matrix before and do not transpose them 

#largest10_var_v_FC_min <- as.data.frame(sort(t(var_v_FC), decreasing = F)[1:10])
largest100_var_v_FC <- sort(var_v_FC, decreasing = T)[1:100]
largest100_var_v_FC= as.matrix(largest100_var_v_FC)

# gives a the 100 genes with highest variance --> may be used as biomarkers 


#################################################################################################

#################################################################################################

### treated data

# creat variance data
var_TreatedVorinostat=(apply(TreatedVorinostat,1,var))
var_TreatedVorinostat

#sort the data

# largest10_var_TreatedVorinostat_min <- sort(var_TreatedVorinostat, decreasing = F)[1:10]
largest100_var_TreatedVorinostat<- sort(var_TreatedVorinostat, decreasing = T)[1:100]
largest100_var_TreatedVorinostat = as.matrix(largest100_var_TreatedVorinostat)

#################################################################################################

#################################################################################################

### compairison of treated and FC results 

# creat vector with gene names 
FC.genes = row.names(largest100_var_v_FC)
Treated.genes = row.names(largest100_var_TreatedVorinostat)

# check if it??s the same length
length(FC.genes)
length(Treated.genes)

setequal(FC.genes,Treated.genes)
# returns FALSE, not the same genes

diff= setdiff(FC.genes,Treated.genes)
length(diff)
# there are 76 genes which are not the same 
# possible reason: houskeeping genes are always upregulated (see treated) 
# and do not change that much by treatment with vorinostat (see FC)
  
