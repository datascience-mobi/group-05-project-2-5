# PCA OF TREATED CELL LINES


# These were my considerations regarding the PCA of treated cell lines,
# as well as the PCA of variables above the 75% quantile 

########################### INSTALLING DATA  #########################    
install.packages("ggfortify")
library(readr)
library(ggplot2)
library("ggfortify") 
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)
t_df_treated = as.data.frame( t(readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))))

########################### PREPARING CELL LINE NAMES AND DRUG NAMES #########################   
# In order to add color the PCA results, we will need a column in our data frames with 
# cell line names OR drug names 
# (see: "https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html" )

# Divide colnames into 1- Cell_line, 2- Drug_name, 3- Dosage, 4- Time_point
m = as.data.frame(strsplit( x= row.names(t_df_treated), split= "_"))
# Make a data frame which consists of a single column with cell_line names
cell_line_names = as.data.frame (t(m[1,]))
# Create a new df from the tranposed df_treated with an additional coloumn of cell_line_names
df_treated_col_cellline = t_df_treated
df_treated_col_cellline$cell_line = cell_line_names[,1]


# Make a data frame which consists of a single column with drug_names
drug_names = as.data.frame (t(m[2,]))
# Create a  data frame from transposed df_treated df with an additional coloumn of drug_names
df_treated_col_drug_name = t_df_treated
df_treated_col_drug_name$drug = drug_names[,1]

########################### PCA OF TREATED CELL LINES #########################    

pca_df_treated = prcomp( df_treated, center =T, scale = T)

autoplot(pca_df_treated)

plot(pca_df_treated , type = "l") # First 2 PC are most relevant

# How to add colours for EACH DRUG 
autoplot(pca_df_treated, data= df_treated_col_drug_name, colour = 'drug')

# How to add colours for EACH CELL LINE
autoplot(pca_df_treated, data= df_treated_col_cellline , colour = 'cell_line')

########################### PCA OF TREATED  > Q75% ########################
# PCA OF TREATED CELL LINES WITH A VARIANCE HIGHER THAN THE 75% OF VARIANCES
 
topVar = apply(t_df_treated, 2, var)

#Select the genes which have a Variance higher than 75% of the variances
df_treated_topVar = t_df_treated[ topVar > quantile(topVar, probs =  0.75)]


df_Tre_topVar_cellline = df_treated_topVar
df_Tre_topVar_cellline$cell_line = cell_line_names[,1]

df_Tre_topVar_drug_name = df_treated_topVar
df_Tre_topVar_drug_name$drug = drug_names[,1]

# Start the PCA  of topVar
pca_df_Tre_topVar = prcomp( df_treated_topVar, center =T, scale = T)
autoplot(pca_df_Tre_topVar)
plot(pca_df_Tre_topVar, type = "l")    # First 5 PCs are most relevant

# How to add colours for EACH DRUG 
autoplot(pca_df_Tre_topVar, data= df_Tre_topVar_drug_name, colour = 'drug')

# How to add colours for EACH CELL LINE
autoplot(pca_df_Tre_topVar, data= df_Tre_topVar_cellline , colour = 'cell_line')

