# PCA OF TREATED CELL LINES

# These were my considerations regarding the PCA of treated cell lines,
# as well as the PCA of variables above the 75% quantile 

install.packages("ggfortify")
library(readr)
library(ggplot2)
library("ggfortify") 
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)
df_treated = as.data.frame( t(readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))))


# In order to add color the PCA results, we will need a column in our data frames with 
# cell line names OR drug names (see: https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html)

# Cut cell_line_names and drug names [l]. 
m = as.data.frame(strsplit( x= row.names(df_treated), split= "_"))
cell_line_names = as.data.frame (t(m[1,]))
drug_names = as.data.frame (t(m[2,]))


# Create 2 different data with an additional coloum each, either with cell_line_names or drug_names
# And order the rows according to cell_line_names or genes_names

df_treated_col_cellline = df_treated
df_treated_col_cellline$cell_line = cell_line_names[,1]
df_treated_col_cellline = df_treated_col_cellline[ order(df_treated_col_cellline$cell_line),]

df_treated_col_drug_name = df_treated
df_treated_col_drug_name$drug = drug_names[,1]
df_treated_col_drug_name = df_treated_col_drug_name[ order(df_treated_col_drug_name$drug),]

# Start the PCA
pca_df_treated = prcomp( df_treated, center =F, scale = F)
# Centering and Scaling are left out, as this could cause 
# an information loss concerning batching effect 
autoplot(pca_df_treated)

plot(pca_df_treated , type = "l") # First 2 PC are most relevant

# How do i get colours for EACH DRUG 
autoplot(pca_df_treated, data= df_treated_col_drug_name, colour = 'drug')

# How do i get colours for EACH CELL LINE
autoplot(pca_df_treated, data= df_treated_col_cellline , colour = 'cell_line')

# PCA OF TREATED CELL LINES WITH A VARIANCE HIGHER THAN THE 75% OF VARIANCES
topVar = apply(df_treated, 2, var)

#Select the genes which have a Variance higher than 75% of the variances
df_treated_topVar = df_treated[ topVar > quantile(topVar, probs =  0.75)]


df_Tre_topVar_col_cellline = df_treated_topVar
df_Tre_topVar_col_cellline$cell_line = cell_line_names[,1]
df_Tre_topVar_col_cellline = df_treated_topVar[ order(df_Tre_topVar_col_cellline$cell_line),]

df_Tre_topVar_col_drug_name = df_treated_topVar
df_Tre_topVar_col_drug_name$drug = drug_names[,1]
df_Tre_topVar_col_drug_name = df_Tre_topVar_col_drug_name[ order(df_Tre_topVar_col_drug_name$drug),]

# Start the PCA  
pca_df_Tre_topVar = prcomp( df_treated_topVar, center =F, scale = F)
autoplot(pca_df_Tre_topVar)
plot(pca_df_Tre_topVar, type = "l")    # First 2 PC are most relevant

# How do i get colours for EACH DRUG 
autoplot(pca_df_Tre_topVar, data= df_Tre_topVar_col_drug_name, colour = 'drug')

# How do i get colours for EACH CELL LINE
autoplot(pca_df_treated, data= df_treated_col_cellline , colour = 'cell_line')

