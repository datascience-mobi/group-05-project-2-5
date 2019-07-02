
#### Exploring correlation in biomarkers using visualization tools #####




### 1.    LOADING DATA                                                                                        ####                                         

## Loading packages
library(readr)
library(rstudioapi)
library(BBmisc)     

## Finding local directory
wd = dirname(rstudioapi::getSourceEditorContext()$path)


## Reading the data


Untreated <- readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))

Treated <- readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))

Sensitivity <- readRDS(paste0(wd,"/data/NegLogGI50.rds"))

Basal <- readRDS(paste0(wd,"/data/CCLE_basalexpression.rds"))

Copynumber <- readRDS(paste0(wd,"/data/CCLE_copynumber.rds"))

Mutations <- readRDS(paste0(wd,"/data/CCLE_mutations.rds"))


Metadata = read.table(paste0(wd,"/data/NCI_TPW_metadata.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)

Cellline_annotation = read.table(paste0(wd,"/data/cellline_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)

Drug_annotation = read.table(paste0(wd,"/data/drug_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)


## Treating the data 
Treated <- as.data.frame(Treated)

Untreated <- as.data.frame(Untreated)

Sensitivity<- as.data.frame(Sensitivity)

## Levels 
levels(Metadata$drug)

##Scaling the data
Treated   = normalize(Treated,   method= "scale")

Untreated = normalize(Untreated, method= "scale")


# Preparing metadata

# Exclude rows, that are treated with other drugs
Metadata_Vorinostat   =  Metadata    [ which(grepl  ("vorinostat" , Metadata$drug) ), ] 
Metadata_Vorinostat   =  Metadata_Vorinostat [ order(Metadata_Vorinostat$tissue), ,drop = FALSE ] 

# => Discover cell line to tissue relations
unique(sort(Metadata$tissue)) # [1] Breast   CNS      Colon    Leukemia Lung     Melanoma Ovarian  Prostate Renal


Metadata_V_after = subset(Metadata_Vorinostat, Metadata_Vorinostat$dose == "5000nM")
Metadata_V_after = Metadata_V_after [, -which(colnames(Metadata_V_after) %in% c( "time", "drug", "sample", "dose" ))]
rownames(Metadata_V_after) = Metadata_V_after$cell
Metadata_V_after
colnames(cor1_tab) = Metadata_V_after$tissue

##  1.1   Creating Vorinostat                                                                                 ####


#Find cell lines, which belong to vorinostat:
#Untreated matrix:
UntreatedVorinostatcolumns <- grep(pattern = "vorinostat",colnames(Untreated))
UntreatedVorinostatcolumns
#-> seems like column 761 - 819 belongs to vorinostat. Check if that is true:
identical(UntreatedVorinostatcolumns, 761:819)
#-> TRUE

#Same with treated matrix:
TreatedVorinostatcolumns <- grep(pattern = "vorinostat",colnames(Treated))
TreatedVorinostatcolumns
identical(TreatedVorinostatcolumns, 761:819)


#Define Vorinostat-data: 
UntreatedVorinostat <- Untreated[,UntreatedVorinostatcolumns]
TreatedVorinostat <- Treated[,TreatedVorinostatcolumns]

#fold change matrix 
FC <- TreatedVorinostat - UntreatedVorinostat

# sensitivity 
vorinostat_Sensitivity_alleZeilen= grep ('vorinostat', rownames(Sensitivity))
vorinostat_Sensitivity_data= Sensitivity[vorinostat_Sensitivity_alleZeilen,]


#We see, that we have a heavy tailed distribution. Nevertheless, we assume normality.




##  1.2   Creating FC Data -  Finding the Biomarkers                                                          ####

FC <- TreatedVorinostat - UntreatedVorinostat

# work with mean of the rows because we only want to compare the genes 
FC_meanrow= rowMeans(FC)


### sort the data 

# work with absolute value to find the highest values
# because we want to have the most up and down regulated genes 
FC_abs= abs(FC_meanrow)

## sort the values to get the 100 largest values 

sortedFC_abs <- sort(FC_abs, decreasing = TRUE)
sortedFC_abs <- as.matrix(sortedFC_abs)

# take the first n for biomarkers 
biomarkers_FC = sortedFC_abs[1:100,]
biomarkers_FC <- as.matrix(biomarkers_FC)


biomarkers_FC30 = sortedFC_abs[1:30,]
biomarkers_FC30 <- as.matrix(biomarkers_FC30)


biomarkers_FC100 = sortedFC_abs[1:100,]
biomarkers_FC100 <- as.matrix(biomarkers_FC100)

### Create matrix with FC values (positive and negativ) 

# add the abs values to FC matrix 
FC_both= cbind(FC_meanrow,FC_abs)
FC_both=as.data.frame(FC_both)

# order this matrix 
FC_both_sorted <- FC_both[order(FC_both$FC_abs, decreasing = TRUE),]

# FC values of biomarkers
# take the first 100 of the sorted matrix, should be the same as un biomarkers_FC 
biomarkers_FC_values = FC_both_sorted[1:100,]

biomarkers_FC_values30 = FC_both_sorted[1:30,]

biomarkers_FC_values100 = FC_both_sorted[1:100,]


# remove the absolute values
biomarkers_FC_values <- subset( biomarkers_FC_values, select = -FC_abs)
biomarkers_FC_values = as.matrix(biomarkers_FC_values)

biomarkers_FC_values30 <- subset( biomarkers_FC_values30, select = -FC_abs)
biomarkers_FC_values30 = as.matrix(biomarkers_FC_values30)

biomarkers_FC_values100 <- subset( biomarkers_FC_values100, select = -FC_abs)
biomarkers_FC_values100 = as.matrix(biomarkers_FC_values100)

# OBJECTIVE DECLARATION


### 2.    CORRELATION 1: Biomarkers found using the mean vs tissue                                            ####

##  2.1   Table with Biomarkers, difference beatween treated and untreated, and cell lines                    ####

## CREATING THE TABLE FOR 30 BIOMARKERS
cor1_tab = FC [ which(row.names(FC) %in% rownames(biomarkers_FC_values30)), ]

#Cleaning the names of the columns and rows
#COLUMNS ONLY WITH CELL LINE + CHANGE NAME OF GENE HIST, CHANGE CELL LINES NAME'S TO TISSUES

colnames(cor1_tab)
col_names_cor1    = as.data.frame(strsplit(x=colnames(cor1_tab),split="_vorinostat"))
colnames (cor1_tab) = as.data.frame (t(col_names_cor1[1,]))[,1]

rownames(cor1_tab)
row_names_cor1    = as.data.frame(strsplit(x=rownames(cor1_tab),split="///HIST"))
rownames (cor1_tab) = as.data.frame (t(row_names_cor1[1,]))[,1]

cor1_tab = as.matrix(cor1_tab)

# This table can only be read if it is a numerical matrix, so first we check if this is true

class(cor1_tab)
# matrix

is.numeric(cor1_tab)
# TRUE

# Now that we have the table, and we know that it is a numeric matrix, we can proceed to use it to
#use for the heatmap

## CREATING THE TABLE FOR 100 BIOMARKERS
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

# This table can only be read if it is a numerical matrix, so first we check if this is true

class(cor1.2_tab)
# matrix

is.numeric(cor1.2_tab)
# TRUE

# Now that we have the table, and we know that it is a numeric matrix, we can proceed to use it to
# use for the heatmap


#   2.2.1 Heatmap with 30 Biomarkers and 59 cell lines                                                        #####

# load packages
library(pheatmap)
library("DESeq")
library(dendextend)


pheatmap(cor1_tab)
cor1.1_tab = cor1_tab

colnames(cor1.1_tab) = Metadata_V_after$tissue

pheatmap(cor1.1_tab)


# Scaling the rows

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

cor1_tab_norm <- t(apply(cor1_tab, 1, cal_z_score))
pheatmap(cor1_tab_norm)


# Ploting a dendogram and cutting the tree
cor1_hclust <- hclust(dist(cor1_tab), method = "complete")

as.dendrogram(cor1_hclust) %>%
  plot(horiz = TRUE)

cor1_hclust_tree <- cutree(tree = as.dendrogram(cor1_hclust), k = 2)

head(cor1_hclust_tree)

# Creating a dataframe, needed for annotations, and adding row annotations

#Cluster annotation
cor1_hclust_tree <- data.frame(cluster = ifelse(test = cor1_hclust_tree == 1, yes = "Cluster 1", no = "Cluster 2"))

head(cor1_hclust_tree)

#Tissue Annotation: Creating a table with the tissue that corresponds to each cell line
tissue = as.data.frame(tissue)
cor1_tissue <- tissue[-c(1:760),] 
cor1_tissue = as.matrix(cor1_tissue)

rownames(cor1_tissue)
row_names_cor1_t    = as.data.frame(strsplit(x=rownames(cor1_tissue),split="_vorinostat"))
colnames(cor1_tissue)[colnames(cor1_tissue)=="X786.0"] <- "786.0"
rownames (cor1_tissue) = as.data.frame (t(row_names_cor1_t[1,]))[,1]

cor1_tissue = as.matrix(t(cor1_tissue))



#Tissue Annotation: Creating the function for the annotation

cor1_tissue <- data.frame(sample = rep(c("Renal", "Lung", "Breast", "Leukemia", "Colon", "Prostate", "Ovarian", "Melanoma", "CNS"), c(4,2)))
row.names(cor1_tissue) <- colnames(cor1_tab)

cor1_tissue <- data.frame(cor1_tissue)
row.names(cor1_tissue) <- colnames(cor1_tab)



meta <- data.frame(
  c(rep("Breast", ncol(cor1_tissue)/9), rep("CNS", ncol(cor1_tissue)/9), rep("Colon", ncol(cor1_tissue)/9), rep("Leukemia", ncol(cor1_tissue)/9),
    rep("Lung", ncol(cor1_tissue)/9), rep("Melanoma", ncol(cor1_tissue)/9), rep("Ovarian", ncol(cor1_tissue)/9), rep("Prostate", ncol(cor1_tissue)/9)
  ),
  row.names=colnames(cor1_tab))
colnames(metadata) <- c("Tissue")




#

cor1_tab_breaks <- seq(min(cor1_tab), max(cor1_tab), length.out = 10)
plot(cor1_tab_breaks)
cor1_tab_breaks

#Heatmap

cor1_colour = list(cluster = c("Cluster 1" = "#68f9f1", "Cluster 2" = "#c9ff87"))


cor1 = pheatmap(cor1_tab,
                annotation_colors = cor1_colour,
                annotation_row = cor1_hclust_tree,
                fontsize = 6.5,
                fontsize_row= 5, 
                fontsize_col = 6,
                gaps_col=50,
                info = TRUE
)



#Heatmap with breaks

pheatmap(cor1_tab,
         annotation_colors = cor1_colour,
         annotation_row = cor1_hclust_tree,
         fontsize = 6.5,
         fontsize_row= 5, 
         fontsize_col = 6,
         gaps_col=50,
         info = TRUE,
         cutree_rows = 2,
         cutree_cols = 2)

pheatmap(cor1_tab,
         annotation_colors = cor1_colour,
         annotation_row = cor1_hclust_tree,
         fontsize = 6.5,
         fontsize_row= 5, 
         fontsize_col = 6,
         gaps_col=50,
         info = TRUE,
         cutree_rows = 2,
         cutree_cols = 4)

# Retrieving hierachical clustering

cor1_ret_hclust <- pheatmap(cor1_tab, silent = TRUE)

cor1_ret_hclust$tree_row %>%
  as.dendrogram() %>%
  plot(horiz = TRUE)


# Cell lines can be observed on the y-axis on the right side, and biomarkers can be observed on the
#x-axis on the bottom. 


#   2.2.2 Heatmap with 100 Biomarkers and 59 cell lines                                                       #####

pheatmap(cor1.2_tab)


# Scaling the rows

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

cor1.2_tab_norm <- t(apply(cor1.2_tab, 1, cal_z_score))
pheatmap(cor1.2_tab_norm)


# Ploting a dendogram and cutting the tree
cor1.2_hclust <- hclust(dist(cor1.2_tab), method = "complete")

as.dendrogram(cor1.2_hclust) %>%
  plot(horiz = TRUE)

cor1.2_hclust_tree <- cutree(tree = as.dendrogram(cor1.2_hclust), k = 2)

head(cor1.2_hclust_tree)

# Creating a dataframe, needed for annotations, and adding row annotations

#Cluster annotation
cor1.2_hclust_tree <- data.frame(cluster = ifelse(test = cor1.2_hclust_tree == 1, yes = "Cluster 1", no = "Cluster 2"))

head(cor1.2_hclust_tree)

#Tissue Annotation: Creating a table with the tissue that corresponds to each cell line
tissue = as.data.frame(tissue)
cor1_tissue <- tissue[-c(1:760),] 
cor1_tissue = as.matrix(cor1_tissue)

rownames(cor1_tissue)
row_names_cor1_t    = as.data.frame(strsplit(x=rownames(cor1_tissue),split="_vorinostat"))
colnames(cor1_tissue)[colnames(cor1_tissue)=="X786.0"] <- "786.0"
rownames (cor1_tissue) = as.data.frame (t(row_names_cor1_t[1,]))[,1]

cor1_tissue = as.matrix(t(cor1_tissue))



#Tissue Annotation: Creating the function for the annotation

cor1_tissue <- data.frame(sample = rep(c("Renal", "Lung", "Breast", "Leukemia", "Colon", "Prostate", "Ovarian", "Melanoma", "CNS"), c(4,2)))
row.names(cor1_tissue) <- colnames(cor1_tab)

cor1_tissue <- data.frame(cor1_tissue)
row.names(cor1_tissue) <- colnames(cor1_tab)



meta <- data.frame(
  c(rep("Breast", ncol(cor1_tissue)/9), rep("CNS", ncol(cor1_tissue)/9), rep("Colon", ncol(cor1_tissue)/9), rep("Leukemia", ncol(cor1_tissue)/9),
    rep("Lung", ncol(cor1_tissue)/9), rep("Melanoma", ncol(cor1_tissue)/9), rep("Ovarian", ncol(cor1_tissue)/9), rep("Prostate", ncol(cor1_tissue)/9)
  ),
  row.names=colnames(cor1_tab))
colnames(metadata) <- c("Tissue")




#

cor1.2_tab_breaks <- seq(min(cor1.2_tab), max(cor1.2_tab), length.out = 10)
plot(cor1.2_tab_breaks)
cor1.2_tab_breaks

#Heatmap for a 100 Biomarkers

cor1.2_colour = list(cluster = c("Cluster 1" = "#68f9f1", "Cluster 2" = "#c9ff87"))


cor1.2 = pheatmap(cor1.2_tab,
                  annotation_colors = cor1.2_colour,
                  annotation_row = cor1.2_hclust_tree,
                  fontsize = 6.5,
                  fontsize_row= 5, 
                  fontsize_col = 6,
                  gaps_col=50,
                  info = TRUE
)



#Heatmap with breaks

pheatmap(cor1.2_tab,
         annotation_colors = cor1.2_colour,
         annotation_row = cor1.2_hclust_tree,
         fontsize = 6.5,
         fontsize_row= 5, 
         fontsize_col = 6,
         gaps_col=50,
         info = TRUE,
         cutree_rows = 2,
         cutree_cols = 3)

pheatmap(cor1.2_tab,
         annotation_colors = cor1.2_colour,
         annotation_row = cor1.2_hclust_tree,
         fontsize = 6.5,
         fontsize_row= 5, 
         fontsize_col = 6,
         gaps_col=50,
         info = TRUE,
         cutree_rows = 2,
         cutree_cols = 5)

# Retrieving hierachical clustering

cor1.2_ret_hclust <- pheatmap(cor1.2_tab, silent = TRUE)

cor1.2_ret_hclust$tree_row %>%
  as.dendrogram() %>%
  plot(horiz = TRUE)


# Cell lines can be observed on the y-axis on the right side, and biomarkers can be observed on the
#x-axis on the bottom. 

##  2.3   Dendogram                                                                                           ####
#   2.3.1 Dendogram with 30 Biomarkers and 59 cell lines                                                      #####

# load packages
library(dendextend)
library(ggplot2)
library(circlize)


## CELL LINES

# ggplot integration

d_gg <- cor1_tab %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=2) %>% set("branches_lwd", c(1.5,1,1)) %>%
  set("branches_lty", c(1,1,3,1,1,2)) %>%
  set("labels_colors") %>% set("labels_cex", c(.5,0.8)) %>% 
  set("nodes_pch", 19) %>% set("nodes_col", c("orange", "black", "plum", NA),
                               main = "How Vorinostat biomarkers define the relatiosnhip between cell lines",
  )

plot(d_gg)

# Mirror horizontal dendogram

hc_hmd <- hclust(dist(cor1_tab), "ave")
d_hmd <- as.dendrogram(hc_hmd)
d_hmd <- d_hmd %>% color_branches(k=2) %>% color_labels
par(mar = c(3,1,1,7))
plot(d_hmd, horiz  = TRUE)

# Polar dendogram

hc_pd <- hclust(dist(cor1_tab))
d_pd <- as.dendrogram(hc_pd)

d_pd <- d_pd %>% 
  color_branches(k=2) %>% 
  color_labels

par(mar = rep(0,4))
circlize_dendrogram(d_pd, labels_track_height = NA, dend_track_height = .3) 

# Retrieving hierachical clustering from the heatmap

cor1_ret_hclust <- pheatmap(cor1_tab, silent = TRUE)

hc_hmd <- hclust(dist(cor1_tab), "ave")
cor1_ret_hclust$tree_row %>% 
  as.dendrogram() %>%
  plot(horiz = TRUE)

cor1_ret_hclust$tree_row = cor1_ret_hclust$tree_row %>% color_branches(k=2) %>% color_labels

par(mar = c(3,1,1,7))
plot(cor1_ret_hclust$tree_row, horiz  = TRUE)


## BIOMARKERS

# ggplot integration

d_gg <- t(cor1_tab) %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=2) %>% set("branches_lwd", c(1.5,1,1)) %>%
  set("branches_lty", c(1,1,3,1,1,2)) %>%
  set("labels_colors") %>% set("labels_cex", c(.5,0.8)) %>% 
  set("nodes_pch", 19) %>% set("nodes_col", c("orange", "black", "plum", NA),
                               main = "How Vorinostat biomarkers define the relatiosnhip between cell lines",
  )

plot(d_gg)

# Mirror horizontal dendogram

hc_hmd <- hclust(dist(t(cor1_tab)), "ave")
d_hmd <- as.dendrogram(hc_hmd)
d_hmd <- d_hmd %>% color_branches(k=2) %>% color_labels
par(mar = c(3,1,1,7))
plot(d_hmd, horiz  = TRUE)

# Polar dendogram

hc_pd <- hclust(dist(t(cor1_tab)))
d_pd <- as.dendrogram(hc_pd)

d_pd <- d_pd %>% 
  color_branches(k=2) %>% 
  color_labels

par(mar = rep(0,4))
circlize_dendrogram(d_pd, labels_track_height = NA, dend_track_height = .3) 


# Retrieving hierachical clustering from the heatmap

cor1_ret_hclust <- pheatmap(t(cor1_tab), silent = TRUE)

hc_hmd <- hclust(dist(t(cor1_tab)), "ave")
cor1_ret_hclust$tree_row %>% 
  as.dendrogram() %>%
  plot(horiz = TRUE)

cor1_ret_hclust$tree_row = cor1_ret_hclust$tree_row %>% color_branches(k=2) %>% color_labels

par(mar = c(3,1,1,7))
plot(cor1_ret_hclust$tree_row, horiz  = TRUE)


#   2.3.2 Dendogram with 100 Biomarkers and 59 cell lines                                                     #####

pheatmap(cor1.2_tab)


# Scaling the rows

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

cor1.2_tab_norm <- t(apply(cor1.2_tab, 1, cal_z_score))
pheatmap(cor1.2_tab_norm)


# Ploting a dendogram and cutting the tree
cor1.2_hclust <- hclust(dist(cor1.2_tab), method = "complete")

as.dendrogram(cor1.2_hclust) %>%
  plot(horiz = TRUE)

cor1.2_hclust_tree <- cutree(tree = as.dendrogram(cor1.2_hclust), k = 2)

head(cor1.2_hclust_tree)



#Heatmap for a 100 Biomarkers

cor1.2_colour = list(cluster = c("Cluster 1" = "#68f9f1", "Cluster 2" = "#c9ff87"))


cor1.2 = pheatmap(cor1.2_tab,
                  annotation_colors = cor1.2_colour,
                  annotation_row = cor1.2_hclust_tree,
                  fontsize = 6.5,
                  fontsize_row= 5, 
                  fontsize_col = 6,
                  gaps_col=50,
                  info = TRUE
)


# Retrieving hierachical clustering

cor1.2_ret_hclust <- pheatmap(cor1.2_tab, silent = TRUE)

cor1.2_ret_hclust$tree_row %>%
  as.dendrogram() %>%
  plot(horiz = TRUE)



#   2.3.3 Comparison of dendograms with 100 Biomarkers and 30 Biomarkers                                      ####
##  2.4   Correlogram and Scatter Plot for 30 Biomarkers and 59 cell lines                                    ####

# Loading Data
library(corrgram)
library(ggplot2)

# Correlogram
corrgram(cor1_tab, order=NULL, panel=panel.shade, text.panel=panel.txt,
         main="Correlogram of Biomarkers and Cell lines") 

# Scatter plot, 3 ways
pairs(cor1_tab, horInd = 1:59, verInd = 1:59, col="#1E90FF")

pairs(cor1_tab, horInd = 1:30, verInd = 1:30, col="#1E90FF")

pairs(cor1_tab, horInd = 1:10, verInd = 1:10, col="#1E90FF")

#Matrix for correlations: with 20 biomarkers and with 100 biomarkers

cor = cor(cor1_tab)
heatmap(cor, col = cm.colors(256))

cor1.2 = cor(cor1.2_tab)
heatmap(cor1.2, col = cm.colors(256))

# Pending 
qplot(cor1_tab, bins = 30)



### 3.    CORRELATION 2: Biomarkers found using the variance vs tissue                                        ####


##2.1 Table with Biomarkers, difference beatwenn treated and untreated, and cell lines
##  3.1   Table with Biomarkers, difference beatween treated and untreated, and cell lines                    ####
##  3.2.1 Heatmap with 30 Biomarkers and 59 cell lines                                                        ####

#   3.2.2 Heatmap with 100 Biomarkers and 59 cell lines                                                       #####
##  3.3   Dendogram                                                                                           ####

#   3.3.1 Dendogram with 30 Biomarkers and 59 cell lines                                                      ####
#   3.3.2 Dendogram with 100 Biomarkers and 59 cell lines                                                     ####
#   3.3.3 Comparison of dendograms with 100 Biomarkers and 30 Biomarkers                                      ####
##  3.4   Correlogram and Scatter Plot for 30 Biomarkers and 59 cell lines                                    ####


### 4.    Comparing the biomarkers obtained using the variance and the mean                                   ####
### 5.    Comparing mean expression of 'biomarker candidates' before and after treatment                      ####

TreatedVorinostat_meanrow= rowMeans(TreatedVorinostat)
UntreatedVorinostat_meanrow= rowMeans(UntreatedVorinostat)


TV_abs= abs(TreatedVorinostat_meanrow)
UV_abs= abs(UntreatedVorinostat_meanrow)

TV_sorted <- sort(TV_abs, decreasing = TRUE)
UV_sorted <- sort(UV_abs, decreasing = TRUE)

TV_matrix <- as.matrix(TV_sorted)
UV_matrix <- as.matrix(UV_sorted)


TV_biomarkers = TV_matrix[1:30,]
UV_biomarkers = UV_matrix[1:30,]

TV_biomarkers <- as.matrix(TV_biomarkers)
UV_biomarkers <- as.matrix(UV_biomarkers)

#Merging the two tables

comp_biomarkers= merge(UV_biomarkers, TV_biomarkers, by = "row.names", all = TRUE)
summary(comp_biomarkers)

# We can observe there are 7 NA's, which means most genes that have a high mean for their 
# expression value before treatment, also have a high mean for their expression value after
# treatment. 
# Therefore, the expression of these genes was already high (high can mean both very upregulated or 
# downregulated) before treatment


# It would be interesting to see if the biomarkers where up- or downregulated after treatment
# Untreated

com_reg_UV = cbind(UntreatedVorinostat_meanrow,UV_sorted)
com_reg_UV = as.data.frame(com_reg_UV)


sort_comp_reg_UV <- com_reg_UV[order(com_reg_UV$UV_sorted, decreasing = TRUE),]


bm_comp_reg_UV = sort_comp_reg_UV[1:30,]


bm_comp_reg_UV <- subset( bm_comp_reg_UV, select = -UV_sorted)
bm_comp_reg_UV = as.matrix(bm_comp_reg_UV)

# Treated

com_reg_TV = cbind(TreatedVorinostat_meanrow,TV_sorted)
com_reg_TV = as.data.frame(com_reg_TV)


sort_comp_reg_TV <- com_reg_TV[order(com_reg_TV$TV_sorted, decreasing = TRUE),]


bm_comp_reg_TV = sort_comp_reg_TV[1:30,]


bm_comp_reg_TV <- subset( bm_comp_reg_TV, select = -TV_sorted)
bm_comp_reg_TV = as.matrix(bm_comp_reg_TV)


#Merge

bm_comp_reg= merge(bm_comp_reg_UV, bm_comp_reg_TV, by = "row.names", all = TRUE)


#Something went wrong and the first ten genes are being displayed, not the idea, I must fix it!!!!

















