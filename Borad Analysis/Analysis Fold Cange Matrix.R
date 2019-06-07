

####################################################################################################

############## Analyse Fold Change Matrix ###############


FC_all = (Treated - Untreated)

############## Mittelwert der Medikamente mit Zelllinien ####################
FC_all_mean=colMeans(FC_all)
FC_all_mean_mat=as.matrix(FC_all_mean)

all_drug_names = Metadata[0:819,]$drug

palette(rainbow(15))
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

###### find Top 10 values 
FC_all_min= apply(FC_all,2,min)
FC_all_max= apply(FC_all,2,max)

largest10_FC_all_min <- sort(FC_all_min, decreasing = F)[1:10]
largest10_FC_all_min = as.matrix(largest10_FC_all_min)
# these genes may be downregulated by the diffrent drugs

largest10_FC_all_max <-sort(FC_all_max, decreasing = T)[1:10]
largest10_FC_all_max = as.matrix(largest10_FC_all_max)
# these genes may be upregulated by the diffrent drugs


#### ideal: can be used to narrow down the values --> which FC values are 'significant' 






