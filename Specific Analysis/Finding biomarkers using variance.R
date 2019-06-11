# FINDING BIOMARKERS USING VARIANCE: (AND COMPARISON THEREOF WITH BIOMARKERS CALCULATED BY T-TEST)

# -> (11.06.2019) I could not work with t-tests from other users, so added my own calculations instead,
# which however are not in their final form

# Complete until the 6th section


# Loading Data                                                                                  ######

install.packages("BBmisc")       #For normalization
library(BBmisc)     

Treated   = normalize(Treated,   method= "scale")
Untreated = normalize(Untreated, method= "scale")

install.packages("dplyr")        #   dplyr for data manipulation
install.packages("ggpubr")       # ggpubr for an easy ggplot2-based data visualization
library("dplyr")
library("ggpubr")
library("dendextend")  


Untreated   = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
Treated     = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
Metadata    = read.table(paste0(wd,"/data/NCI_TPW_metadata.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)


# Drug sensitivity assay
Sensitivity = readRDS(paste0(wd,"/data/NegLogGI50.rds"))

# Basal molecular profiles of cancer cell lines 
Basal       = readRDS(paste0(wd,"/data/CCLE_basalexpression.rds"))
Copynumber  = readRDS(paste0(wd,"/data/CCLE_copynumber.rds"))
Mutations   = readRDS(paste0(wd,"/data/CCLE_mutations.rds"))

# Feature annotation
# cell line metadata => To be expanded with Basal Molecular Profiles
Cellline_annotation = read.table(paste0(wd,"/data/cellline_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)
# mechanism of action, targets ect.
Drug_annotation = read.table(paste0(wd,"/data/drug_annotation.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)



# As Vorinostat is our chosen drug, we should cut the coulumn names including "vorinostat"
Vorinostat_Untreated   =  Untreated [, which(grepl  ( "vorinostat" ,colnames(Untreated) )  )] 
Vorinostat_Treated     =  Treated   [, which(grepl  ( "vorinostat" ,colnames(Treated) )  )] 

# Calculate the amount of the difference in expression levels => Drug Response
Drug_Response =  Vorinostat_Treated - Vorinostat_Untreated

Vorinostat_Untreated

# Cleaning the column names
col_names_VU    = as.data.frame(strsplit(x=colnames(Vorinostat_Untreated),split="_"))
colnames (Vorinostat_Untreated) = as.data.frame (t(col_names_VU[1,]))[,1]


col_names_VT    = as.data.frame(strsplit(x=colnames(Vorinostat_Treated),split="_"))
colnames (Vorinostat_Treated) = as.data.frame (t(col_names_VT[1,]))[,1]




# 0-  I will compare the variances of each gene before and after treatment among 59 cell lines. ####

# Genes showing high variance before treatment are either cancerogenic
# or have high variance due to cell line specificity

# Thus, I will detect the genes with highly variant activity among 59 cell lines
# and measure the difference in variances before and after treatmen (Variance_after - Variance_before)

# I expect the genes with highly variant activity (before treatment), which also show a high variance difference after 
# treatment to be responding to the drug

# For genes, who show a high variance difference (Variance_after - Variance_before), 
# but did not show high variance before treatment: 
# - These could be either omnipresent cancerogenic genes with less variance before treatment,
#   since they play roles in general mechanisms such as proliferation or 
#   apoptosis supression, and reacted differently to the drug among cell lines.
#   ( However, this does not seem very likely to me )
# - Or side-effects on non-cancerogenic genes

# In 4, I will look at the variance of gene responses to drug.
# Genes expression levels change differently for a gene among cell lines.
# A high variance in drug response would indicate that for a cancerogenic gene in a specific 
# cell line, its aberrant expression levels normalized, causing different amount
# of change and thus variance. 



# 1-  Genes with highest variance before treatment                                              ####
Variances_Untreated = as.data.frame(  apply(Vorinostat_Untreated,1,var)  )
colnames(Variances_Untreated)[1] = "var"

Most_Variances_Unt = Variances_Untreated[ order(-Variances_Untreated$var), , drop= FALSE ] 

nintyfive_quartile_Untr = quantile(Variances_Untreated$var, probs = 0.95)

Most_Variances_Unt = Most_Variances_Unt[-which(Most_Variances_Unt$var < nintyfive_quartile_Untr),, drop= FALSE] 
dim(Most_Variances_Unt)    #  [1] 665   1

#  These are the genes showing highest variance among untreated cell lines
#  and can be considered as relevant for the cancerogenic activity

# 2-  Variance of genes after treatment                                                         ####

# Later, I will compare variances before and adter treatment and see if there is a change. 
# A change would signify a drug response.

Variances_Treated = as.data.frame(  apply(Vorinostat_Treated,1,var)  )
colnames(Variances_Treated)[1] = "var"

Most_Variances_Tr = Variances_Treated[ order(-Variances_Treated$var), ,drop = FALSE ] 
colnames(Most_Variances_Tr)[1] = "var"

nintyfive_quartile_Tr = quantile(Variances_Treated$var, probs = 0.95)

Most_Variances_Tr = Most_Variances_Tr[-which(Most_Variances_Tr$var < nintyfive_quartile_Tr),, drop= FALSE] 
dim(Most_Variances_Tr)    #  [1] 665   1

# 665 genes from both data are not the same
rownames(Most_Variances_Unt) %in% rownames(Most_Variances_Tr)

# 3-  Variance Comparison                                                                       ####
# 3a- High Variance change in genes who showed high variance before treatment                   ####
Variance_Change = as.data.frame (Variances_Treated[,1] - Variances_Untreated[,1])

row.names(Variance_Change)   = rownames(Variances_Treated)
colnames(Variance_Change)[1] = "Var_Change"
Variance_Change = Variance_Change[ order(-Variance_Change$Var_Change), ,drop = FALSE ] 

upper_quartile_VC   = quantile(Variance_Change$Var_Change, probs = 0.975)
lower_quartile_VC   = quantile(Variance_Change$Var_Change, probs = 0.025)

Most_pos_var_change  = Variance_Change[-which(Variance_Change$Var_Change < upper_quartile_VC),, drop= FALSE] 
Most_neg_var_change  = Variance_Change[-which(Variance_Change$Var_Change > lower_quartile_VC),, drop= FALSE] 
Most_neg_var_change  = Most_neg_var_change[ order(Most_neg_var_change$Var_Change), ,drop = FALSE ] 

Highest_VC = append(rownames(Most_neg_var_change) , rownames(Most_neg_var_change))

# Which genes whose variance changed strongly were in the highly variant group before treatment
Biomarker_candidates = intersect( Highest_VC, rownames(Most_Variances_Unt)  )
Biomarker_candidates # 242 out of 666


# 3b- High Variance change in genes whose variance was not very high in the beginning           ####

# These are genes whose variance changed strongly, but were not variant before treatment.
# These may be cancerogenic genes (equally expressed in all cell lines) or housekeeping genes,
# whose expression level was differently affected by drug in different cell lines.

# In each case, we cannot come to a relevant conclusion.


# 4-  Variance in Drug Response                                                                 ####

# After the treatment, we expect different levels of drug response for a gene 
# across different cell lines if the gene is responsible for an oncological activity.

# If there is less variance in drug response, we could say that the gene is 
# equally affected by drug in different cell lines and is therefore not specific for a cancer type.
Variances_DR = as.data.frame(  apply(Drug_Response,1,var)  )
colnames(Variances_DR)[1] = "var"

Variances_DR = Variances_DR[ order(-Variances_DR$var), , drop= FALSE ] 
nintyfive_quartile_VC  = quantile(Variances_DR$var, probs = 0.95)

Most_Variances_DR = Variances_DR[-which(Variances_DR$var < nintyfive_quartile_VC),, drop= FALSE ]
hist(Most_Variances_DR[,1])

# 5-  Comparison of genes from 3 and 4    => Biomarkers (calculated by variance)                ####

# In 3, I selected the genes, which I think are responsible for oncological character
# and also respond to the treatment.

# In 4, I focused on the genes with the highest variance in drug response.
# Higher variance in drug response means, for specific cancer types a gene was overexpressed
# only in some cell lines (therefore high variance before treatment), and the response is
# also variant among cell lines, as some cell lines are affected differently than others.
Biomarker_candidates     #     247 
nrow(Most_Variances_DR)  # [1] 665
Biomarkers= as.data.frame(intersect( Biomarker_candidates, rownames(Most_Variances_DR)) ) 
Biomarkers  # 171 biomarkers, this number will decrease after the comparison with t-test values !!!!!

rownames(Biomarkers)     =  Biomarkers$`intersect(Biomarker_candidates, rownames(Most_Variances_DR))`
# This long name was name of thecolumn =>   $`intersect(Biomarker_candidates, rownames(Most_Variances_DR))`


# Now I create data frame Vorinostat_Treated and Vorinostat_Untreated only with the biomarker genes
Biomarkers_Untreated     =  Vorinostat_Untreated[ which(row.names(Vorinostat_Untreated) %in% rownames(Biomarkers)),]
Biomarkers_Treated       =  Vorinostat_Treated  [ which(row.names(Vorinostat_Treated) %in% rownames(Biomarkers)),  ]
Drug_Response_Biomarkers =  Biomarkers_Treated - Biomarkers_Untreated
# And also a Drug Response file with only biomarker genes
Drug_Response_Biomarkers


# 6-  Calculating Biomarkers using t-test                                                       ####

# We will calculate the genes which show smallest p-values:
# -first,  using a paired two-tailed t-test on treated and basal expression levels
# -second, using a independent two-tailed t-test on drug responses (gene expression change levels)

# We will then compare the genes showing drug response with very high significance (highest p values)


# 6a- Calculating Biomarkers using t-test ( paired; Treated vs Untreated)                       ####
p_values_Vorinostat = sapply(1:13299, function(i){  
  t.test( Vorinostat_Treated[i,],Vorinostat_Untreated[i,])$p.value
})

p_values_Vorinostat           = as.data.frame (p_values_Vorinostat)
rownames(p_values_Vorinostat) = rownames((Vorinostat_Treated))
colnames(p_values_Vorinostat) =  " p_values"

p_values_Vorinostat


Significance_Treated = sapply(1:13299, function(i) {
  # Return whether p-value  is < 0.05 
  Vorinostat_p_Values[i,] < 0.025
})

p_values_Vorinostat$significance = Significance_Treated
p_values_Vorinostat = p_values_Vorinostat[ which (p_values_Vorinostat$significance == TRUE),]
p_values_Vorinostat[1:5,]
p_values_Vorinostat = p_values_Vorinostat[ order(p_values_Vorinostat$` p_values`), ,drop = FALSE ]

hist(p_values_Vorinostat[,1])

p_values_Vorinostat = head( p_values_Vorinostat,30, drop =FALSE   )
p_values_Vorinostat
p_values_Vorinostat[1:10,]



# 6b- Calculating Biomarkers using t-test ( independent; Drug Respomse = Treated - Untreated)   ####

Expr_change_Vorinostat =  Vorinostat_Treated - Vorinostat_Untreated
Expr_change_Vorinostat = normalize(Expr_change_Vorinostat, method= "scale")



t.test(Expr_change_Vorinostat[1,])$p.value

p_values_EC = sapply(1:13299, function(i){  
  t.test( Expr_change_Vorinostat[i,])$p.value
})

p_values_EC            = as.data.frame (p_values_EC )
rownames(p_values_EC ) = rownames((Vorinostat_Treated))
colnames(p_values_EC ) =  " p_values"


Significance_EC = sapply(1:13299, function(i) {
  # Return whether p-value  is < 0.05 
  p_values_EC[i,] < 0.05
})
Significance_EC

p_values_EC$significance = Significance_EC
p_values_EC = p_values_EC[ which (p_values_EC$significance == TRUE),]
dim(p_values_EC)
p_values_EC = p_values_EC[ order(p_values_EC$` p_values`), ,drop = FALSE ]
p_values_EC
hist(p_values_EC[,1])


p_values_EC =  head(p_values_EC, 30, drop = FALSE)
p_values_EC[1:50,]






# 7-  Comparing Biomarkers calculated by variance and those calculated by t-test               ####
# 8-  PCA of biomarkers        = Dissecting Patterns of Transcriptome Regulation                ####

Biomarkers_Untreated     =  Vorinostat_Untreated[ which(row.names(Vorinostat_Untreated) %in% rownames(Biomarkers)),  ]
Biomarkers_Treated       =  Vorinostat_Treated  [ which(row.names(Vorinostat_Treated) %in% rownames(Biomarkers)),  ]
Drug_Response_Biomarkers = Biomarkers_Treated - Biomarkers_Untreated


pca = prcomp(Drug_Response_Biomarkers, center = T, scale. = T)  

plot(pca, type = "l")
# Make sure to have library(ggfortify)
autoplot(pca)


# Now we can group cell lines using omics files and color the cell lines belonging to the same group 

# 9-  kmeans of biomarkers     = Dissecting Patterns of Transcriptome Regulation                ####

t_DR_Biom = as.data.frame( t(Drug_Response_Biomarkers)  )
t_DR_Biom

topVar = apply(Drug_Response_Biomarkers, 1, var)
summary(topVar)

topVar_DR_Biom = Drug_Response_Biomarkers[topVar > quantile(topVar, probs = 0.75), ]
dim(topVar_DR_Biom)

rm(topVar)


# Creating a correlation based distance matrix
cor.mat = cor(topVar_DR_Biom, method = "spearman")
cor.dist = as.dist(1 - cor.mat)
cor.hc = hclust(cor.dist, method = "ward.D2")
cor.hc = as.dendrogram(cor.hc)

plot(cor.hc, las = 2, cex.lab = 0.7) 


# 10- Working with Anootations and Metadata                                                     ####

# Experiment metadata for the gene expression profiling above. 
Metadata
# => Discover cell line to tissue relations
unique(sort(Metadata$tissue)) # [1] Breast   CNS      Colon    Leukemia Lung     Melanoma Ovarian  Prostate Renal

Metadata_BM   =  Metadata    [ which(grepl  ("vorinostat" , Metadata$drug) ), ] 
Metadata_BM   =  Metadata_BM [ order(Metadata_BM$tissue), ,drop = FALSE ] 

sort(Metadata_BM$cell)  # For every cell line, there is a treated and 
Metadata_BM$dose
Metadata_BM   =  Metadata_BM [ -which(Metadata_BM$dose == "0nM"), ] 

Metadata_BM

Metadata_BM$cell %in% Biomarker_candidates



# Now we can colour our PCA results from 6

autoplot(pca)
autoplot(pca, data = Metadata, colour = 'tissue')


# We cannot, because here we have only tissue to cell line connections, not the genes

sapply ( )



# Drug sensitivity assay = NetLogGI50
Sensitivity  

# Basal molecular profiles of cancer cell lines 
Basal       
Copynumber 
Mutations   

# Feature annotation
# cell line metadata => To be expanded with Basal Molecular Profiles
Cellline_annotation 
# mechanism of action, targets ect.
Drug_annotation 










