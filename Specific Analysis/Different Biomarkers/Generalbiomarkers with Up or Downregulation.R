#You first need to load the "General biomarkers(absolute values)" file!!!

#We want to add the information, whether the FC of genes is mainly positive or negative, to the generalbiomarker-matrix:
#Matrix with information about general up- or downregulation:
FCVorinostatmean_noabs_matrix <- as.matrix(FCVorinostatmean_noabs)
#(In this matrix, FCmean was calculated not with absolute FC values)

#Vector with genenames of the top100 biomarkers:
generalbiomarkergenes


#loop, which adds those biomarkergenes to a vector, which are mostly upregulated:
i=1
Upregulated <- c()
while(i<=100) {
  if(FCVorinostatmean_noabs_matrix[generalbiomarkergenes[i],1] >0) {
    Upregulated <- union(Upregulated, generalbiomarkergenes[i])
  }
  i=i+1}

#Create a vector that includes for each gene sequentially if it is up- or downregulated:
i=1
Generalchange <- c()
while(i<=100) {
  if(FCVorinostatmean_noabs_matrix[generalbiomarkergenes[i],1] >0) {
    Generalchange <- append(Generalchange, "Up")
  } else {
    Generalchange <- append(Generalchange, "Down")
  }
  i=i+1}

#Bind the information about up/downregulation as a new column to the biomarkermatrix:
top100generalbiomarkers_withUporDown <- cbind(top100generalbiomarkers, Generalchange)
colnames(top100generalbiomarkers_withUporDown)[1] <- "FCmean"

#--> the new biomarkermatrix: top100generalbiomarkers_withUporDown


