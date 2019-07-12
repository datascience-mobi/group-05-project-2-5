# comparison genes with most variance and biomarkers 

VarianceTreated <- apply(Treated_norm, 1, var)
VarianceTreatedSorted <- sort(VarianceTreated, decreasing = TRUE)
VarianceTreatedSorted=as.data.frame(VarianceTreatedSorted[1:100])
a=rownames(VarianceTreatedSorted)

generalbiomarkergenes=as.vector(generalbiomarkergenes)
setequal(a, generalbiomarkergenes)
diff=setdiff(a, generalbiomarkergenes)
length(diff)

#--> genes with most variance in whole data not in our biomarkers 
