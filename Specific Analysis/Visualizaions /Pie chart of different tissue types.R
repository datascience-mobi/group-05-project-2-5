library(readr)
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)

#Load data:

Treated = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
Metadata = read_tsv(paste0(wd,"/data/NCI_TPW_metadata.tsv"))

#Use tissue information from Metadata and add it to Treated-matrix, then select those columns which belong to Vorinostattreatment:

Metadata <- as.matrix(Metadata)
Metadatatissue <- Metadata[1:819,"tissue"]
Treatedwithtissue <- rbind(Treated, Metadatatissue)

    TreatedVorinostatcolumns <- grep(pattern = "vorinostat",colnames(Treatedwithtissue))
    TreatedVorinostattissue <- Treatedwithtissue[,TreatedVorinostatcolumns]

#Create a pie chart:
    
palette(rainbow(9))
pie(table(TreatedVorinostattissue["Metadatatissue",]),
    main="Tissue types of Vorinostat-treated samples", 
    radius=1, 
    cex=0.8)

# better coloring: http://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization

