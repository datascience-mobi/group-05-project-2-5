# PCA OF TREATED CELL LINES

# 0-  INSTALLING DATA                                                                #########################    
install.packages("ggfortify")
library(readr)
library(ggplot2)
library("ggfortify") 
library(rstudioapi)

wd = dirname(rstudioapi::getSourceEditorContext()$path)
Metadata = read.table(paste0(wd,"/data/NCI_TPW_metadata.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
Treated = as.data.frame(readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds")))


# 1-  PREPARING CELL LINE NAMES AND DRUG NAMES                                       #########################

# In order to add color the PCA results, we will need a column in our data frames with
# tissue names OR drug names
# (see: "https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html" )

# 1a- Extracting tissue names from Metada                                            #########

# Metadata and Treated include partly the same experiments
colnames(Treated) == Metadata$sample

# Reduce metadata to experiment included in "Treated" data
Metadata   =  as.data.frame (Metadata [which(  Metadata$sample == colnames(Treated)  ),])

# Transpose the Treated data
t_Treated = as.data.frame( t(Treated) )
dim(t_Treated)

# Add a new column with Metadata names
t_Treated$tissue = Metadata$tissue
head(Metadata$tissue,30)
dim(t_Treated)
head(t_Treated[,13300])

# Now, we have at the end of each row (demonstrating a cell line) the corresponding tissue name


# 1b- Extracting drug names from experiment names                                    #########

# Surely, this part could also be done with ( t_Treated$tissue = Metadata$tissue) function
# But I wanted to keep the way I used at the beginning for nostalgic reasons.


# Divide colnames into 1- Cell_line, 2- Drug_name, 3- Dosage, 4- Time_point
m = as.data.frame(strsplit( x= row.names(t_Treated), split= "_"))
m[1:4,1:4]
# Make a data frame which consists of a single column with drug_names
drug_names = as.data.frame (t(m[2,]))
head(drug_names,62)
# Add a column to the data frame t_Treated an additional coloumn of drug_names
dim(t_Treated)
t_Treated$drug_names = drug_names[,1]

# Did it work according to the plan?
tail(colnames( t_Treated  ),4)     # Yes, last columns are:   [1] "ZZEF1"   "ZZZ3"   "tissue"  "drug_names"

# Does each column contain what they are supposed to? 
head (t_Treated$tissue, 7)         # Yes,  [1] Renal    Lung     Breast   Renal    Leukemia Prostate Colon   
head (t_Treated$drug_names,3)      #       [1] 5-Azacytidine 5-Azacytidine 5-Azacytidine

t_Treated[,13300:13301] 




# 2a- PCA OF TREATED,  ALL CELL LINES                                                #########################    
dim(t_Treated)
# PCA of each gene (column), except the last 2, which contain tissue names and drug names
pca_Treated = prcomp( t_Treated[,1:13299], center =F, scale = F)

plot(pca_Treated , type = "l")      # First 2 PC are most relevant#
autoplot(pca_Treated)               # There are different clusters, 

# so now we should look if these clusters are relevant to drug or tissue:


# HOW TO ADD COLOURS FOR EACH DRUG 
autoplot(pca_Treated, data= t_Treated, colour = 'drug_names')

# We can see that different drugs cause different clusters:
# Now we can check if this has to do with batch effect !

pca_normalized_Treated = prcomp( t_Treated[,1:13299], center =T, scale = T)
autoplot(pca_normalized_Treated, data= t_Treated, colour = 'drug_names')

# Scaling and centering clearly changes the shape of our data
# Also the importance (in percent) of PCs change very strongly !

plot(pca_normalized_Treated, type = "l")  # First 3-5 PCs could be relevant



# HOW TO ADD COLURS FOR EACH CELL LINE
autoplot(pca_Treated, data= t_Treated , colour = 'tissue')   # Looks chaotic

autoplot(pca_normalized_Treated, data= t_Treated, colour = 'tissue') 
# There are interrupting structures 


# FOCUSING ON THE MOST VARIANT STRUCTURE, WE CAN SEE THE EFFECTS OF DRUG PERTURBANCE MORE CLEARLY.



# 2b- PCA OF TREATED,  VARIANCE  >  Q75%                                             ####

# PCA OF TREATED CELL LINES WITH A VARIANCE HIGHER THAN THE 75% OF VARIANCES
dim(t_Treated)

# Calculating variances of genes
Var_Treated = as.data.frame (apply(t_Treated[,1:13299], 2, var))
dim(Var_Treated)

# What is the 75% quantile
third_quantile = quantile(Var_Treated[,1], probs =  0.75)
third_quantile

# Select the genes which have a Variance higher than 75% of the variance
topVar_Treated = as.data.frame(t_Treated[, Var_Treated > third_quantile], drop = FALSE)
dim(topVar_Treated)       # [1]  819 3326

# Is 3326 one fourth of 13299 ?
13299/4      # [1] 3324.75
# So we should have 3325 + 2 columns ( "tissue" and "drug_name")

# What is missing?
topVar_Treated[1:3,3323:3326]     #  ZSCAN31     ZWINT      ZYX    drug_names
                                  #  No tissue column !

"tissue" %in% colnames(topVar_Treated) # [1] FALSE

# Add the tissue column:
topVar_Treated$tissue = Metadata$tissue

topVar_Treated[1:3,3323:3327]     # Now it works



# Start the PCA  of topVar
pca_topVar_Treated = prcomp( topVar_Treated[,1:3325], center =T, scale = T )
autoplot(pca_topVar_Treated)
plot(pca_topVar_Treated, type = "l")    # First 5 PCs are most relevant



# How to add colours for EACH DRUG 
autoplot(pca_topVar_Treated, data= topVar_Treated, colour = 'drug_names')

# Comparing different PCs
autoplot(pca_topVar_Treated, x=1, y=3, data= topVar_Treated, colour = 'drug_names')
autoplot(pca_topVar_Treated, x=1, y=4, data= topVar_Treated, colour = 'drug_names')
autoplot(pca_topVar_Treated, x=1, y=5, data= topVar_Treated, colour = 'drug_names')

autoplot(pca_topVar_Treated, x=2, y=3, data= topVar_Treated, colour = 'drug_names')
autoplot(pca_topVar_Treated, x=2, y=4, data= topVar_Treated, colour = 'drug_names')
autoplot(pca_topVar_Treated, x=2, y=5, data= topVar_Treated, colour = 'drug_names')

autoplot(pca_topVar_Treated, x=3, y=4, data= topVar_Treated, colour = 'drug_names')
autoplot(pca_topVar_Treated, x=3, y=5, data= topVar_Treated, colour = 'drug_names')

autoplot(pca_topVar_Treated, x=4, y=5, data= topVar_Treated, colour = 'drug_names')


# How to add colours for EACH TISSUE
autoplot(pca_df_Tre_topVar, data= topVar_Treated , colour = 'tissue')

# Comparing different PCs
autoplot(pca_topVar_Treated, x=1, y=3, data= topVar_Treated , colour = 'tissue')
# Leukemia and lung are distinguishable ( PC1 (10.06%), PC3 (6.94%) )
autoplot(pca_topVar_Treated, x=1, y=4, data= topVar_Treated , colour = 'tissue')
autoplot(pca_topVar_Treated, x=1, y=5, data= topVar_Treated , colour = 'tissue')

autoplot(pca_topVar_Treated, x=2, y=3, data= topVar_Treated , colour = 'tissue')
autoplot(pca_topVar_Treated, x=2, y=4, data= topVar_Treated , colour = 'tissue')
autoplot(pca_topVar_Treated, x=2, y=5, data= topVar_Treated , colour = 'tissue')


autoplot(pca_topVar_Treated, x=3, y=4, data= topVar_Treated , colour = 'tissue')
# Leukemia, renal and melanoma can be distinguished( PC3 (6.94%), PC4 (5.44%)  )
autoplot(pca_topVar_Treated, x=3, y=5, data= topVar_Treated , colour = 'tissue')
# Leukemia, renal and melanoma can be distinguished( PC3 (6.94%), PC5 (3.5%)  )


autoplot(pca_topVar_Treated, x=4, y=5, data= topVar_Treated , colour = 'tissue')
# Leukemia is a clearly different group in this comparison ( PC4 (5.44%), PC5 (3.5%)  )
