
# copy biomarker and Biomarker3 from annotation to data 
# barplot with diffrent tissues 

library(readr)
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)


vorinostat_annotation=read.csv2(paste0(wd,"/../data/biomarker.csv"),header = TRUE, quote="\"")

tissue=vorinostat_annotation$affected.Tissue..if.specific.

col=palette(rainbow(6))
table(tissue)
barplot(table(tissue), ylab="counts", main="affected tissues by biomarker", col=col)


######### with nones 
vorinostat_annotation=read.csv2(paste0(wd,"/../data/Biomarker3.csv"),header = TRUE, quote="\"")

tissue=vorinostat_annotation$affected.Tissue..if.specific.

col=palette(rainbow(7))
table(tissue)
barplot(table(tissue), ylab="counts", main="affected tissues by biomarker", col=col)