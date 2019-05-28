####################################################################################################

############## Analyse Fold Change Matrix ###############


FC_all = (Treated - Untreated)

############## Mittelwert der Medikamente mit Zelllinien ####################
FC_all_mean=colMeans(FC_all)
FC_all_mean_mat=as.matrix(FC_all_mean)

all_drug_names = Metadata[0:819,]$drug

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
df_test <- data.frame("ID" = 1:819, "Color" = color_array, "Name" = all_drug_names, "MeanValue" = as.vector(FC_all_mean_mat))

#
barplot( height = df_test$MeanValue
         , names.arg = df_test$Name
         , las = 3
         , col = df_test$Color
         , border = NA
)


############ scatter plot ##############
plot(FC_all_mean, col= Metadata$tissue)
plot(FC_all_mean, col= Metadata$drug)


########## denity plot ###############
plot(density(FC_all_mean))



####################################################################################################

####################################################################################################

######### Einfuegen der Medikamentenzeile #############

#Treated_colnames = Treated[0,]

# Auftrennen der Spaltenbezeichnungen 
gesplittete_col_names = as.data.frame(strsplit(x=colnames(Treated),split="_"))

# Laenge auf treated "normieren"
names(Treated) = names(gesplittete_col_names[2,])

# verbinden 
Treated_with_names2= rbind(Treated ,gesplittete_col_names[2,]) 





