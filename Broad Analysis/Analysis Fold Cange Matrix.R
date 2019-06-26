
################################# Analyse Fold Change Matrix #####################################

##################################################################################################

##################################################################################################

########################################### old way ##############################################


# creat FC data 
FC_all = (Treated - Untreated)
FC_all_mean=colMeans(FC_all)
#FC_all_mean_mat=as.matrix(FC_all_mean)

# creat names vector 
all_drug_names = Metadata$drug
length(all_drug_names) 

# prepare barplot 

#farben die spaeter im plot benutzt werden
color_vector_rainbow = rainbow(15)
#index fuer das colorvector_rainbow array
index = -1
#string um zu checken, ob eine neue categorie ist
old_name = ""
#array indem spaeter die farben gespeichert werden
color_array = c()
#loop ueber alle namen
for(i in 0:length(all_drug_names)){
  #print(paste(all_drug_names[i], " - ",old_name, "Bool:", !identical(old_name, all_drug_names[i])))
  # Wenn der gen name sich aendert zaehlen wir den index hoch, da wir dann eine andere farben haben moechten
  if(!identical(old_name, all_drug_names[i])){
    index = index + 1
    old_name = all_drug_names[i]
  }
  # fuer jeden namen fuege die aktuelle farbe hinzu
  color_array[i] <- color_vector_rainbow[index]
}
# erstelle dataframe mit den spalten 
df_test <- data.frame("ID" = 1:819, "Color" = color_array, "Name" = all_drug_names, "MeanValue" = as.vector(FC_all_mean))

#barplot 
barplot( height = df_test$MeanValue, names= FALSE, col = df_test$Color, border = NA)

##################################################################################################

##################################################################################################

########################################### new way ##############################################

# creat FC data 
FC_all = (Treated - Untreated)
FC_all_mean=colMeans(FC_all)

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









