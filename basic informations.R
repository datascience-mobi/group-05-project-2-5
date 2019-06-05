################################################################################################
####GSVA 
#installation
https://bioconductor.org/packages/release/bioc/html/GSVA.html
#paper mit den plots (siehe auch Vortrag)
https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf 
#code zum paper
https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.R
#manual
https://bioconductor.org/packages/release/bioc/manuals/GSVA/man/GSVA.pdf

#### Enhanced volcano --> besser dokumentiert als GSVA volcano 
https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

#### infos zu KEGG
https://bioconductor.statistik.tu-dortmund.de/packages/3.1/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.pdf
# hab ich selbst noch nicht angeschaut 
(((https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)))


################################################################################################

################ volcano plot ##################################################################

# erster versuch, Problem: beide Packages nicht installierbar 

# Moeglichkeit 1
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

# Moeglichkeit 2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")

# enhanced volcano braucht matrix mit 'normalen' daten und p.values 
pValues <- apply(VorinostatTotal, 1, function(x) t.test(x[col_untreated],x[col_treated],paired = TRUE, alternative = "two.sided")$p.value)
pValues.data = as.data.frame(pValues)

test=cbind(FC,pValues.data)
#vielleicht kann das programm nur mit zwei spalten arbeiten 
test2=cbind(FC_meanrow,pValues)

View(test)
EnhancedVolcano(test,
                +                 lab = rownames(test),
                +                 x = 'log2FoldChange',
                +                 y = 'pvalue',
                +                 xlim = c(-5, 8))