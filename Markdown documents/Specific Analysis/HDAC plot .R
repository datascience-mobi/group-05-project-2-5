# search for HDAC in biomarkers 
target <- grep(pattern = "HDAC",colnames(top100generalbiomarkers))
length(target)

target_vorinostat=c("HDAC1","HDAC10","HDAC11","HDAC2","HDAC3","HDAC5","HDAC6",
                    "HDAC8","HDAC9")

FC_meanrow = as.data.frame(FC_meanrow)
genes = row.names(FC_meanrow)
FC_new=cbind(FC_meanrow,genes)



# filter data 
library(dplyr)
FC_target = as.data.frame(filter(FC_new, genes %in% target_vorinostat))
FC_target

##### boxplots 

# vector with  biomarker and target 
last10biomarkers= generalbiomarkergenes[1:10]

c=union(target_vorinostat,last10biomarkers)
length(c)

# find values in FC matrix 
FC= TreatedVorinostat-UntreatedVorinostat
FC=as.data.frame(FC)

FC_HDAC <-FC[c,]

# not all HDAC targets can be found in our data, remove the missing ones 
FC_HDAC =na.omit(FC_HDAC)


# boxplot 
col=palette(rainbow(17))
boxplot(t(FC_HDAC), las=3, col=col, ylab="Fold Change Values", main="Compairison of Biomarker and target gene expression change ")
abline(v=7.5, col="blue", lty=5)


