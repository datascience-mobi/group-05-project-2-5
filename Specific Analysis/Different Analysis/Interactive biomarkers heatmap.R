### Interactive biomarkers heatmap
  
### 1. Load packages                                                                            ####
library(Seurat)
library(hexbin)
library(devtools)
library(heatmaply)


### 2. Load Biomarkers data                                                                     ####

  
#Load creat vorinostat and biomarkers over FC
  
#Matrix for 100 biomarkers
cor1.2_tab = FC [ which(row.names(FC) %in% rownames(biomarkers_FC_values100)), ]
  
#Cleaning the names of the columns and rows
#COLUMNS ONLY WITH CELL LINE + CHANGE NAME OF GENE HIST, CHANGE CELL LINES NAME'S TO TISSUES
  
colnames(cor1.2_tab)
col_names_cor1.2    = as.data.frame(strsplit(x=colnames(cor1.2_tab),split="_vorinostat"))
colnames (cor1.2_tab) = as.data.frame (t(col_names_cor1.2[1,]))[,1]
  
rownames(cor1.2_tab)
row_names_cor1.2    = as.data.frame(strsplit(x=rownames(cor1.2_tab),split="///HIST"))
rownames (cor1.2_tab) = as.data.frame (t(row_names_cor1.2[1,]))[,1]
  
cor1.2_tab = as.matrix(cor1.2_tab)
  
  
#Heatmap 
heatmaply(cor1.2_tab, k_row = 2, k_col = 2)
  
  
### 3. Heatmap                                                                                  #### 
heatmaply(cor1.2_tab, k_row = 2, k_col = 2)


heatmaply(cor1.2_tab, 
          column_text_angle = 90,
          k_row = 3, 
          k_col = 3, 
          main = "Interactive heatmap of 100 biomarkers and 59 cancer cell lines",
          fontsize_row = 6,
          fontsize_col = 6,
          )


### 4. Correlation matrix                                                                       ####

heatmaply(cor(cor1.2_tab), 
          column_text_angle = 90,
          k_row = 3, 
          k_col = 3, 
          main = "Correlation matrox showing 59 cancer cell lines",
          fontsize_row = 6,
          fontsize_col = 6,
          )
