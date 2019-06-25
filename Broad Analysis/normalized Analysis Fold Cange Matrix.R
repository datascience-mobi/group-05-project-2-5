
################################# Analyse Fold Change Matrix #####################################

# creat FC data 
FC_all = (Treated - Untreated)
#normalization
FC_all_norm <- apply(FC_all, 2, function(x){
  (x - mean(x)) / sd(x)
})
FC_all_mean=colMeans(FC_all_norm)

# creat levels for coloring 
drug <- Metadata$drug
palette(rainbow(15))

# creat boxplot
barplot( height = FC_all_mean, names= FALSE, col = drug, border = NA,main= "Fold change")

# creat legend 
levels <- as.factor(levels(drug))
legend("topright", inset = c(-0.3,0), legend= levels(drug), xpd = TRUE, pch=19, col = levels, title = "drugs")

##################################################################################################

##################################################################################################

######################################### scatter plot ###########################################

### drug
plot(FC_all_mean, col= Metadata$drug, main="Fold change with drugs")
# legend 
drug <- Metadata$drug
levels <- as.factor(levels(drug))
legend("topright", inset = c(-0.3,0), legend= levels(drug), xpd = TRUE, pch=19, col = levels, title = "drugs")
# in general similar FC values, only 5-Azacytidine and bortezomib have clear outliners 

### tissue 
plot(FC_all_mean, col= Metadata$tissue,main="Fold change with tissues")
# legend 
tissue <- Metadata$tissue
levels <- as.factor(levels(tissue))
legend("topright", inset = c(-0.3,0), legend= levels(tissue), xpd = TRUE, pch=19, col = levels, title = "tissues")
# no correlation with tissue 




##################################################################################################

##################################################################################################

######################################### density plot ###########################################

plot(density(FC_all_mean), main= "density Fold change")

# normally distributed 
# most values between -0.75 and 0.5 

##################################################################################################

##################################################################################################


#################################### find Top 10 values  #########################################

### find all min and max 
FC_all_min= (apply(FC_all,2,min))
FC_all_max= (apply(FC_all,2,max))

### sort min and max 

# most down regulated gene
largest10_FC_all_min <- (sort(FC_all_min, decreasing = F)[1:10])
largest10_FC_all_min 
# frequent drugs: vorinostat (3), bortezomib (3)
# fequent cell line: OVCAR-4 (3)

# most up regulated genes 
largest10_FC_all_max <- (sort(FC_all_max, decreasing = T)[1:10])
largest10_FC_all_max
# frequent drugs: bortezomib (6)
# fequent cell line: none









