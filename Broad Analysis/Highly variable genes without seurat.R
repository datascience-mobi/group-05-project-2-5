

VarianceTreated <- apply(Treated_norm, 1, var)
VarianceTreatedSorted <- sort(VarianceTreated, decreasing = TRUE)
VarianceTreatedSorted[1:10]

VarianceUntreated <- apply(Untreated_norm, 1, var)
VarianceUntreatedSorted <- sort(VarianceUntreated, decreasing = TRUE)
VarianceUntreatedSorted[1:10]

###TREATED PLOT

x <- c(1:13299)
VarianceTreated <- as.data.frame(VarianceTreated)
VarianceTreatedwithnumbers <- cbind(VarianceTreated, x)
VarianceTreatedwithnumbers <- as.data.frame(VarianceTreatedwithnumbers)


subsettreated <- subset(VarianceTreatedwithnumbers, 
       VarianceTreatedwithnumbers$VarianceTreated >= VarianceTreatedSorted[10])

overlap <- intersect(names(VarianceTreatedSorted[1:10]), names(VarianceUntreatedSorted[1:10]))

subsetboth_treated <- VarianceTreatedwithnumbers[overlap,]




p <- ggplot(VarianceTreatedwithnumbers, aes (x= x, y= VarianceTreated ))

p+ geom_point() +ggtitle("Genes with highest variance in treated data")+
  geom_text(data=subsettreated, label= rownames(subsettreated),hjust=1, vjust=1, size =3)+
  geom_point(data = subsettreated, colour="red", size = 3)+
  geom_point(data= subsetboth_treated, size = 3, color = "orange")  
                      

### UNTREATED PLOT

VarianceUntreated <- as.data.frame(VarianceUntreated)
VarianceUntreatedwithnumbers <- cbind(VarianceUntreated, x)
VarianceUntreatedwithnumbers <- as.data.frame(VarianceUntreatedwithnumbers)


subsetuntreated <- subset(VarianceUntreatedwithnumbers, 
                 VarianceUntreatedwithnumbers$VarianceUntreated >= VarianceUntreatedSorted[10]) 

subsetboth_untreated <- VarianceUntreatedwithnumbers[overlap,]

p <- ggplot(VarianceUntreatedwithnumbers, aes (x= x, y= VarianceUntreated ))

p+ geom_point() +ggtitle("Genes with highest variance in untreated data")+
  geom_text(data= subsetuntreated, label= rownames(subsetuntreated),hjust=0, vjust=0, size =3)+
  geom_point(data = subsetuntreated, colour="red", size = 3)+
  geom_point(data= subsetboth_untreated, size = 3, color = "yellow")

                      

          