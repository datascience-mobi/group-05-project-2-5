# barplot general changes 
# load general biomarkers 

top100generalbiomarkers_withUporDown=as.data.frame(top100generalbiomarkers_withUporDown)
change=top100generalbiomarkers_withUporDown$Generalchange
col=palette(rainbow(2))
barplot(table(a), ylab="counts", main="number of up and down regulated genes", col=col)