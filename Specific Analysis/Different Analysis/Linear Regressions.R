

### Linear Regressions ###

## Change log?

## Normalization? 

##### PART 1 ##################################################################################################################
####  1.  LOADING DATA                                                                                                      ####
###   1.1 Packages                                                                                                          ####
###   1.2 Biomarkers                                                                                                        ####
##  (1)  Creating Vorinostat


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




##  (2)  Creating FC Data -  Finding the Biomarkers

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


####  2.  PREPARING THE DATA                                                                                                ####
###   2.1 log change                                                                                                        ####
###   2.2 Creating the necessary tables (Copynumber with biomarkers/all genes and Cell lines with drug sensitivity)         ####

## (1) Table 1: Selection of 100 Biomarkers in copynumber 

BM_copynumber = Copynumber[ which(row.names(Copynumber) %in% rownames(biomarkers_FC_values100)), ]

BM_Copynumber_meancol= colMeans(Copynumber)

CN = as.data.frame(BM_Copynumber_meancol)


## (2) Table 2: All genes in copynumber 

CN_meancol= colMeans(Copynumber)

CN_all = as.data.frame(CN_meancol)


## (3) Table 3: Selection of Doubling time from cellline_annotation

library(dplyr)

# Selecting the desired columns

Doubling_Time <- cellline_annotation %>%
  select(Cell_Line_Name, Doubling_Time)

# Changig the name of the "name" column to the names of the cell lines

row.names(Doubling_Time) <- Doubling_Time$Cell_Line_Name
Doubling_Time[1] <- NULL

Doubling_Time


####  3.  SIMPLE LINEAR REGRESSION WITH 100 BIOMARKERS: Drug sensitivity with doubling time                                 ####

# Can we predict drug sensitivity using doubling time?
# How much of the variance of the data can be explained using the doubling time?

###   3.1 Linear Regression                                                                                                 ####

# l.g = lm(predicted ~ predictor, data = dat) 

DT = as.data.frame(Doubling_Time)

DS = as.data.frame(drug_sensitivity)

lm_tab = transform(merge(DT,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

reg1 <- lm(vorinostat ~ Doubling_Time, data = lm_tab)

summary(reg1)

###   3.2 Checking the normalization of residuals                                                                           ####

hist(reg1$residuals, breaks = 20)

#data does NOT look normally distributed

qqnorm(reg1$residuals)
qqline(reg1$residuals)

###   3.3 Visualization: Plots                                                                                              ####

plot(lm_tab$vorinostat, lm_tab$Doubling_Time)

#abline(, col = "red", lwd = 2)

####  4.  SIMPLE LINEAR REGRESSION WITH 100 BIOMARKERS: Drug sensitivity with copynumber                                    ####

# Can we predict drug sensitivity using the copynumber data?
# How much of the variance of the data can be explained using the copynumber data?

###   4.1 Linear Regression                                                                                                 ####

CN = as.data.frame(BM_Copynumber_meancol)

DS = as.data.frame(drug_sensitivity)

lm_tab2 = transform(merge(CN,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

reg2 <- lm(vorinostat ~ BM_Copynumber_meancol, data = lm_tab2)

summary(reg2)

###   4.2 Checking the normalization of residuals                                                                           ####

hist(reg2$residuals, breaks = 20)

#data does NOT look normally distributed

qqnorm(reg2$residuals)
qqline(reg2$residuals)


###   4.3 Visualization: Plots                                                                                              #### 
####  5.  SIMPLE LINEAR REGRESSION WITH ALL GENES: Drug sensitivity with copynumber                                         ####

# Can we predict drug sensitivity using the copynumber data?
# How much of the variance of the data can be explained using the copynumber data?

###   5.1 Linear Regression                                                                                                 ####
                                                                                                                          #### 
CN_all = as.data.frame(CN_meancol)

DS = as.data.frame(drug_sensitivity)

lm_tab_all2 = transform(merge(CN_all,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

reg_all2 <- lm(vorinostat ~ CN_meancol, data = lm_tab_all2)

summary(reg_all2)

###   5.2 Checking the normalization of residuals                                                                           ####

hist(reg_all2$residuals, breaks = 20)

#data does NOT look normally distributed

qqnorm(reg2$residuals)
qqline(reg2$residuals)
###   5.3 Visualization: Plots                                                                                              ####
####  6.  MULTIPLE LINEAR REGRESSION WITH 100 BIOMARKERS: Drug sensitivity with doubling time and copynumber                ####


###   6.1 Linear Regression                                                                                                 ####

CN = as.data.frame(BM_Copynumber_meancol)

DT = as.data.frame(Doubling_Time)

DS = as.data.frame(drug_sensitivity)

lm_tab_m = transform(merge(CN, lm_tab,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_m <- na.omit(lm_tab_m)

reg_m <- lm(vorinostat ~ BM_Copynumber_meancol + Doubling_Time, data = lm_tab_m)

summary(reg_m)
###   6.2 Checking the normalization of residuals                                                                           ####
###   6.3 Visualization: Plots                                                                                              ####

plot(lm_tab_m)

####  7.  General Conclusions                                                                                               ####










##### PART 2 ##################################################################################################################
####  8.  Using PCA to determine independent variables                                                                      ####
###   8.1 PCA                                                                                                               ####

CN = as.data.frame(BM_Copynumber_meancol)

DT = as.data.frame(Doubling_Time)

DS = as.data.frame(drug_sensitivity)

lm_tab_m = transform(merge(CN, lm_tab,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_m <- na.omit(lm_tab_m)

reg_m <- lm(vorinostat ~ BM_Copynumber_meancol + Doubling_Time, data = lm_tab_m)

summary(reg_m)

lm_tab_pca = lm_tab_m[,c(3,1,2)]

pca = prcomp(lm_tab_pca[, -1])
summary(pca)

par(las = 2)
par(mar = c(1,11,1,2))
barplot(pca$rotation[, 1], horiz = TRUE, main = "PC1", col = "red")



###   8.2 Checking for correlation                                                                                          ####
###   8.3 Linear Regression                                                                                                 ####
###   8.4 Visualization: Plots                                                                                              ####




