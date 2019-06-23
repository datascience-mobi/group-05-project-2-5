

### Interactive heatmap

#Load packages

library(Seurat)
library(hexbin)
library(devtools)
library(heatmaply)

#Heatmap 
heatmaply(cor1.2_tab, k_row = 2, k_col = 2)


heatmaply(cor1.2_tab, 
          column_text_angle = 90,
          k_row = 3, 
          k_col = 3, 
          main = "Interactive heatmap of 100 biomarkers and 59 cancer cell lines",
          fontsize_row = 6,
          fontsize_col = 6,
          )


# Correlation matrix

heatmaply(cor(cor1.2_tab), 
          column_text_angle = 90,
          k_row = 3, 
          k_col = 3, 
          main = "Interactive heatmap of 100 biomarkers and 59 cancer cell lines",
          fontsize_row = 6,
          fontsize_col = 6,
          )

