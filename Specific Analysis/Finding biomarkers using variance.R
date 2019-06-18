# FINDING BIOMARKERS 

# Loading Data                                                                                    ######

install.packages("BBmisc")       #For normalization
library(BBmisc)     

install.packages("dplyr")        #   dplyr for data manipulation
install.packages("ggpubr")       # ggpubr for an easy ggplot2-based data visualization
library("dplyr")
library("ggpubr")
library("dendextend")  


Untreated   = readRDS(paste0(wd,"/data/NCI_TPW_gep_untreated.rds"))
Treated     = readRDS(paste0(wd,"/data/NCI_TPW_gep_treated.rds"))
Metadata    = read.table(paste0(wd,"/data/NCI_TPW_metadata.tsv"), header = TRUE, sep ="\t", stringsAsFactors = TRUE)

Untreated = as.data.frame(Untreated)
Treated = as.data.frame(Treated)

Treated   = normalize(Treated,   method= "scale")
Untreated = normalize(Untreated, method= "scale")

colnames(Treated)
colnames(Untreated)

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

?normalize

# As Vorinostat is our chosen drug, we should cut the coulumn names including "vorinostat"
Vorinostat_Untreated   =  as.data.frame (Untreated [, which(grepl  ( "vorinostat" ,colnames(Untreated) )  )] )
Vorinostat_Treated     =  Treated   [, which(grepl  ( "vorinostat" ,colnames(Treated) )  )] 

# Calculate the amount of the difference in expression levels => Drug Response
Drug_Response =  Vorinostat_Treated - Vorinostat_Untreated

colnames(Vorinostat_Untreated)

# Cleaning the column names
colnames(Vorinostat_Untreated)
col_names_VU    = as.data.frame(strsplit(x=colnames(Vorinostat_Untreated),split="_"))
colnames (Vorinostat_Untreated) = as.data.frame (t(col_names_VU[1,]))[,1]
Vorinostat_Untreated[1:5,1:5]

col_names_VT    = as.data.frame(strsplit(x=colnames(Vorinostat_Treated),split="_"))
colnames (Vorinostat_Treated) = as.data.frame (t(col_names_VT[1,]))[,1]



# 0-      Different type of Biomarkers                                                            ####
# 0.1 -   Biomarkers as genes highly expressed in all cancerogenic cell lines                     ####

# PART 1 

# Here, we will look at the drug responses (Gene expr. after Treatment - Gene expr. before Tr.)
# And see which genes showed the highest responses in average


# 0.2 -   Biomarkers as genes highly expressed in some specific cancerogenic cell lines           ####

# PART 2


# For this part, we will look at variances 
# Firstly, to the genes which show a high variance before treatment => Cancer type specific genes
# Not only cancerogenic genes will show high variance but also cell line specific genes
# (For example, in lung cells vs intestinal cells)

# Then find the high changes in variance => If the variance change, the genes would be responding
# We expect the highest responses to come from cancerogenic genes


#  HOW TO FIND THESE GENES:

# I will detect the genes with highly variant activity among 59 cell lines
# and measure the difference in variances before and after treatmen (Variance_after - Variance_before)

# I expect the genes with highly variant activity (before treatment), 
# which also show a high variance difference after treatment to be responding to the drug


# For genes, who show a high variance difference (Variance_after - Variance_before), 
# but did not show high variance before treatment: 
# - These could be either omnipresent cancerogenic genes with less variance before treatment,
#   since they play roles in general mechanisms such as proliferation or 
#   apoptosis supression, and reacted differently to the drug among cell lines.
#   ( These could be investigated by t-tests )
# - Or side-effects on non-cancerogenic genes


# In part 4, I will look at the variance of gene responses to drug.
# Genes expression levels change differently for a gene among cell lines.
# A high variance in drug response would indicate that for a cancerogenic gene in a specific 
# cell line, its aberrant expression levels normalized, causing different amount
# of change and thus variance. 


# 0.3 -   Why we cannot depend on t-test for finding biomarkers                                   ####     
# 1 -     PART 1    =>         *Biomarkers_Highest_Mean_Drug_Response*                            ####

# (This is the adviced way of finding biomarker by our tutor)

# In this part, we will look at the highest mean values of changes 
# in gene expression levels after treatment
# => Genes whose expression level changed strongly in many cell types

# Mean values of drug response
mean_Drug_Response = as.data.frame(apply(Drug_Response,1, mean))

# We want the highest changes, regardless of their positivity or negativity
# => We will use absolute values
mean_Drug_Response = abs(mean_Drug_Response)

# Now list in decreasing order
mean_Drug_Response = mean_Drug_Response[ order(-mean_Drug_Response[,1]), , drop= FALSE ]

# Get the hundert highes changes
Biomarkers_Highest_Mean_Drug_Response = mean_Drug_Response[ order(-mean_Drug_Response)[1:100], , drop = FALSE ]



# 2 -     PART 2    =>         *Biomarkers_Variance*                                              ####

# The genes found in the first part are those, which showed the highest change on average 
# among 59 cell lines.

# Yet, in our cell lines we have different cancer types and in order to find genes specific
# for different cancer types, we need to look at the genes with high variance.

# The problem will be that not only genes with cancerogenic activity will differ among cell lines,
# but also the tissue specific genes will have high variance (compare: lung versus intestine)

# In order to find the *cancerogenic* genes with high variance, we will
# 1- Look at all the genes with high variance before treatment
# 2- See which of those responded strongly to the tratment
# 3- Expect the strongest responses will be made by cancerogenic genes, as these
# tend to be affected most strongly by chemotherapeutic treatment.

# VARIANCE IN DRUG RESPONSE => CELL LINE SPECIFITY
# Then, we will investigate the highest variances in Drug Response,
# the higher variance in drug response, the more specific the drug response to cell line


# 2.1 -   VARIANCE AFTER - VARIANCE BEFORE => VARIANCE CHANGE                                     ####
# 2.1.1-  Genes with highest variance before treatment                                            ####
Variances_Untreated = as.data.frame(  apply(Vorinostat_Untreated,1,var)  )
colnames(Variances_Untreated)[1] = "var"

Most_Variances_Unt = Variances_Untreated[ order(-Variances_Untreated$var), , drop= FALSE ] 

nintyfive_quartile_Untr = quantile(Variances_Untreated$var, probs = 0.95)

Most_Variances_Unt = Most_Variances_Unt[-which(Most_Variances_Unt$var < nintyfive_quartile_Untr),, drop= FALSE] 
dim(Most_Variances_Unt)    #  [1] 665   1

#  These are the genes showing highest variance among untreated cell lines
#  and can be considered as relevant for the cancerogenic activity

# 2.1.2-  Variance of genes after treatment                                                       ####

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

# 2.1.3-  Variance Comparison   ↓ ↓ ↓                                                             ####
# 2.1.3a- High Variance change in genes which showed high variance before treatment               ####
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


# 2.1.3b- High Variance change in genes which showed LOW variance before treatment                ####

# These are genes whose variance changed strongly, but were not variant before treatment.
# These may be cancerogenic genes (equally expressed in all cell lines) or housekeeping genes,
# whose expression level was differently affected by drug in different cell lines.

# In each case, we cannot come to a relevant conclusion.


# 2.2-    VARIANCE IN DRUG RESPONSE                                                               ####

# After the treatment, we expect different levels of drug response for a gene 
# across different cell lines if the gene is responsible for an oncological activity.

# If there is less variance in drug response, we could say that the gene is 
# equally affected by drug in different cell lines and is therefore not specific for a cancer type.
Variances_DR = as.data.frame(  apply(Drug_Response,1,var)  )
colnames(Variances_DR)[1] = "var"

Variances_DR = Variances_DR[ order(-Variances_DR$var), , drop= FALSE ] 
nintyfive_quartile_VC  = quantile(Variances_DR$var, probs = 0.95)

Most_Variances_DR = Variances_DR[-which(Variances_DR$var < nintyfive_quartile_VC),, drop= FALSE ]
rownames(Most_Variances_DR)

# 2.3-    Comparison of Biomarkers from 3 and 4 =>  Biomarkers_Variance*                          ####

# In 3, I selected the genes, which I think are responsible for oncological character
# and also respond to the treatment.

# In 4, I focused on the genes with the highest variance in drug response.

# Higher variance in drug response means, that a gene responded differently to the treatment
# in different cell lines. This would be a sign for differing drug responses. 

# These genes could be cancerogenic or not, however I expect the cancerogenic genes to show 
# higher variance, as regular genes probably would not have reason to be affected differently
# by chemotherapic agents. 
# (Chemotherapic agents tend to attack the highly dividing cells the most strongly)

Biomarker_candidates     #     247 
nrow(Most_Variances_DR)  # [1] 665
Biomarkers_1= as.data.frame(intersect( Biomarker_candidates, rownames(Most_Variances_DR)) ) 
Biomarkers_1  # 169 biomarkers, this number will decrease after the comparison with t-test values !

# Changing row names
rownames(Biomarkers)     =  Biomarkers$`intersect(Biomarker_candidates, rownames(Most_Variances_DR))`
# This long name was name of the column name:
# =>   $`intersect(Biomarker_candidates, rownames(Most_Variances_DR))`


# Now I create data frame Vorinostat_Treated and Vorinostat_Untreated only with the biomarker genes
Biomarkers_Untreated     =  Vorinostat_Untreated[ which(row.names(Vorinostat_Untreated) %in% rownames(Biomarkers)),]
Biomarkers_Treated       =  Vorinostat_Treated  [ which(row.names(Vorinostat_Treated) %in% rownames(Biomarkers)),  ]
Biomarkers_Variance      =  Biomarkers_Treated - Biomarkers_Untreated

dim(Biomarkers_Variance)   #  [1] 169  59


# 3 -     PART 3     =>        NO Biomarkers_t_test data, see: 3.2                                ####

# t-TEST
# We will also look at the p values of changes, however as the changes are really small,
# a lot of tiny changes will be considered significant, although they are not necessarily so.

# 3.1-    How to calculate Biomarkers using t-test and why                                        ####

# We will calculate the genes showing drug response with very high significance (lowest p values):

# -first,  using a paired two-tailed t-test on treated and basal expression levels
# -second, using an independent two-tailed t-test on drug responses (gene expression level changes)
# lastly, we will compare genes from both t-tests



# 3.1a-   Calculating Biomarkers using t-test (two-tailed paired;      Treated vs. Untreated)     ####

p_values_Vorinostat = sapply(1:13299, function(i){  
 t.test( Vorinostat_Treated[i,],Vorinostat_Untreated[i,])$p.value
})

p_values_Vorinostat           = as.data.frame (p_values_Vorinostat)
rownames(p_values_Vorinostat) = rownames((Vorinostat_Treated))
colnames(p_values_Vorinostat) =  " p_values"

sort(p_values_Vorinostat[,1])  #  A lot of values are essentially zero

hist(p_values_Vorinostat[,1])


Significance_Treated = sapply(1:13299, function(i) {
  # Return whether p-value  is < 0.05  
  Vorinostat_p_Values[i,] < 0.05
})

p_values_Vorinostat$significance = Significance_Treated
p_values_Vorinostat = p_values_Vorinostat[ which (p_values_Vorinostat$significance == TRUE),]
p_values_Vorinostat[1:5,]
p_values_Vorinostat = p_values_Vorinostat[ order(p_values_Vorinostat$` p_values`), ,drop = FALSE ]

hist(p_values_Vorinostat[,1])


Biomarker_candidates_t_test_1 = p_values_Vorinostat

dim(Biomarker_candidates_t_test_1)    #   [1] 665   2



# 3.1b-   Calculating Biomarkers using t-test (two-tailed independent; Treated - Untreated)       ####

Drug_Response = Vorinostat_Treated - Vorinostat_Untreated
Drug_Response = normalize(Drug_Response, method= "scale")


p_values_EC = sapply(1:13299, function(i){  
  t.test(Drug_Response[i,])$p.value
})

p_values_EC            = as.data.frame (p_values_EC )
rownames(p_values_EC ) = rownames((Vorinostat_Treated))
colnames(p_values_EC ) = " p_values"

dim(p_values_EC)

hist(p_values_EC[,1])

# Which values are significant?
Significance_EC = sapply(1:13299, function(i) {
  # Return whether p-value  is < 0.05
p_values_EC[i,] < 0.05
})
Significance_EC

# Create a column with significance = TRUE/FALSE
p_values_EC$significance = Significance_EC
# Delete unsignificant rows
p_values_EC = p_values_EC[ which (p_values_EC$significance == TRUE),]
dim(p_values_EC)   #   [1] 8480    2

# Order data frame according to p values
p_values_EC = p_values_EC[ order(p_values_EC$` p_values`), ,drop = FALSE ]
p_values_EC
hist(p_values_EC[,1])

Biomarker_candidates_t_test_2 =  p_values_EC

dim (Biomarker_candidates_t_test_2)  # [1] 8480    2 => TOO MUCH



# 3.2-    Comparing Biomarkers calculated using t-test                                            ####

# Which of the different biomarker candiates from both t tests are common ?

Biomarkers_t_test  = intersect(rownames(Biomarker_candidates_t_test_1), rownames(Biomarker_candidates_t_test_2))
Biomarkers_t_test  # 665 names

# CONCLUSION: AS EXPECTED, WHEN VALUES ARE TOO SMALL (AS IS THE CASE IN DRUG RESPONSES), 
# SLIGHTEST CHANGES ARE CONSIDERED SIGNIFICANT. (SEE: 3.1b)
# THEREFORE, WE CANNOT WORK WITH T-TEST RESULTS, AS THESE ARE BIASED !






####

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










