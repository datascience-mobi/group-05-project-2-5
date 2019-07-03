
# copy Biomarkers from annotation to data 
# barplot with diffrent tissues 

library(readr)
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)


vorinostat_annotation=read.csv2(paste0(wd,"/../data/Biomarkers.csv"),header = TRUE, quote="\"")

tissue=vorinostat_annotation$affected.Tissue..if.specific.



col=palette(rainbow(9))
table(tissue)
barplot(table(tissue), ylab="counts", main="affected tissues by biomarker", col=col, las=3)


######### with nones 
vorinostat_annotation=read.csv2(paste0(wd,"/../data/Biomarker3.csv"),header = TRUE, quote="\"")

tissue=vorinostat_annotation$affected.Tissue..if.specific.

col=palette(rainbow(7))
table(tissue)
barplot(table(tissue), ylab="counts", main="affected tissues by biomarker", col=col, las=3)

#### barplot functions 
general.function=vorinostat_annotation$general.Function

col=palette(rainbow(20))
table(general.function) 
barplot(table(general.function), ylab="counts", main="general function biomarker", col=col, las=3)


