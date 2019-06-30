

### Linear Regressions ###

## Pending: Change log?
## No linear relationship: lowess, GAM
## Complete until 6

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

BM_Copynumber_meancol= colMeans(BM_copynumber)

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

## Data frames
DT = as.data.frame(Doubling_Time)

DS = as.data.frame(drug_sensitivity)

## Table with drug sensitivity and doubling time per cell line

lm_tab = transform(merge(DT,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)


###   3.1 Plots and visualization: Predicting how fit linear regression will be as a model to describe our data             #### 


## Checking for correlation 

cor(lm_tab)

## Ploting the data: can a linear relationship be observed? Should we expect a high value for R-squared?

# Visualization of doubling time vs drug sensitivity

# (1) Scatter Plot

scatter.smooth(lm_tab$vorinostat, 
               lm_tab$Doubling_Time, 
               col = "dodgerblue1",
               main = "Drug sensitivity & Doubling time Regression",
               xlab = "Drug sensitivity",
               ylab = "Doubling time",
               cex = 1.3,
               pch = 1)

#At first look, the points in the plot are so scattered, that a linear relationship seems unlikely. 
#The second problem that one can observe in this graphic, is that the line describing the behaviour
#is not straight. 

# (2) Box plot
par(mfrow=c(1, 2))
boxplot(lm_tab$vorinostat, main="Drug sensitivity", sub=paste("Outlier rows: ", boxplot.stats(lm_tab$vorinostat)$out)) 
boxplot(lm_tab$Doubling_Time, main="Doubling time", sub=paste("Outlier rows: ", boxplot.stats(lm_tab$Doubling_Time)$out)) 

#A boxplot can help us visualize the amount of outliers in our data. This is relevant as too many (extreme) outliers can
#have great impact on the results of our analysis and can change the outcome completely. They can easily affect the slope. 
#Some outliers can be observed.

# (3) Density: Should be expect normality for drug sensitivity?

library(e1071)

par(mfrow=c(1, 2)) 

plot(density(lm_tab$vorinostat), 
     main="Density Plot: Drug Sensitivity", 
     ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(lm_tab$vorinostat), 2))
     )

polygon(density(lm_tab$vorinostat), col="royalblue1")

#Skewness: 0.22 -> Plot is slightly skewed to the right

plot(density(lm_tab$Doubling_Time), 
     main="Density Plot: Doubling time", 
     ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(lm_tab$Doubling_Time), 2))
    )

polygon(density(lm_tab$Doubling_Time), col="skyblue1")

#Skewness: 1.03 -> Plot is skewed to the right. 

# (4) Correlation: what is the level level of linear dependence between the two variables?

cor(lm_tab$vorinostat, lm_tab$Doubling_Time)

#[1] -0.2057795
#A good value for correlation lies close to 1 or -1, whilst the value 0 is undesirable. Values closer to 0 indicate
#that there is a weak relationship.  


#These analysis help us predcit whether a linear regression is or not the best model to describe our data. 
#Taking into consideration all results so far for this part, it is not unreasonable to predict that a linear
#regression will probably not be the best model to describe the relationships in our data.

###   3.2 Linear Regression                                                                                                 ####

# l.g = lm(predicted ~ predictor, data = dat) 


## Linear Regression

reg1 <- lm(vorinostat ~ Doubling_Time, data = lm_tab)

## Details about the linear regression: what we need draw some conclusions

summary(reg1)

# Multiple R-squared:  0.04235 -> This indicates that only 4,235% percent of the variation in the data (drug sensitivity)
#can be explained by the relationship between drug sensitivity and doubling time. In other words, there is a 4.235% 
#variance reduction when we take the doubling time into account. 

# p-value: 0.1116 
#As the p-value for reg1 is significantly larger than 0.05 and R-squared tells us the doubling time only explains 4.235% 
#of the variation in the data, it is safe to assume that there is no linear relationship between drug sensitivity and 
#doubling time, a.k.a doubling time cannot predict drug sensitivity. 

# More information about the fit (linear ecuation: y = y-intercept + slope * x) : 
confint(reg1)


###   3.3 Checking the normalization of residuals                                                                           ####

hist(reg1$residuals, 
     breaks = 20, 
     xlab = "Residuals", 
     main = "Drug sensitivity vs doubling time: Histogram of the residuals")

# The data does NOT look normally distributed: As the residuals are not normally distributed, then the hypothesis that
#they are a random dataset, takes the value NO.
#This means that the linear regression model does not explain all trends in the dataset.

qqnorm(reg1$residuals)
qqline(reg1$residuals, col = "red")


###   3.4 Visualization: Plots that describe the linear regression                                                          ####

## Residual diagnostics: are the various assumptions that underpin linear regression reasonable for our data?

library(lattice)

xyplot(resid(reg1) ~ fitted(reg1),
       xlab = "Fitted Values",
       ylab = "Residuals",
       main = "Residual Diagnostic Plot",
       col = "slateblue3",
       panel = function(x, y, ...)
       {
         panel.grid(h = -1, v = -1)
         panel.abline(h = 0)
         panel.xyplot(x, y, ...)
       }
      )


qqmath( ~ resid(reg1),
      xlab = "Theoretical Quantiles",
      ylab = "Residuals",
      abline(a=-0.01034, b=1))

#Here we expect to be able to draw a straight line that is fitting for most points. Although, this graphic seems promising,
#too many points would not be properly fit. 

# Visualization of regression 

par(mar = c(4, 4, 2, 2), mfrow = c(1, 2))
plot(reg1, which = c(1, 2))

####  4.  SIMPLE LINEAR REGRESSION WITH 100 BIOMARKERS: Drug sensitivity with copynumber                                    ####

# Can we predict drug sensitivity using the copynumber data?
# How much of the variance of the data can be explained using the copynumber data?

## Data frames

CN = as.data.frame(BM_Copynumber_meancol)

DS = as.data.frame(drug_sensitivity)

## Table with drug sensitivity and doubling time per cell line

lm_tab2 = transform(merge(CN,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
lm_tab2 <- na.omit(lm_tab2)


###   4.1 Plots and visualization: Predicting how fit linear regression will be as a model to describe our data             #### 

## Checking for correlation

cor(lm_tab2)

## Ploting the data: can a linear relationship be observed? Should we have expected a high value for R-squared?

# Visualization of doubling time vs drug sensitivity

# (1) Scatter Plot
scatter.smooth(lm_tab2$vorinostat, 
              lm_tab2$BM_Copynumber_meancol, 
              col = "yellowgreen",
              main = "Drug sensitivity & Copynumber Regression",
              xlab = "Drug sensitivity",
              ylab = "Copynumber",
              cex = 1.3,
              pch = 1)

#At first look, the points in the plot are so scattered, that a linear relationship seems unlikely.
#The second problem that one can observe in this graphic, is that the line describing the behaviour
#is not straight. 


# (2) Box plot

par(mfrow=c(1, 2))
boxplot(lm_tab2$vorinostat, main="Drug sensitivity", sub=paste("Outlier rows: ", boxplot.stats(lm_tab2$vorinostat)$out)) 
boxplot(lm_tab2$BM_Copynumber_meancol, main="Copynumber", sub=paste("Outlier rows: ", boxplot.stats(lm_tab2$BM_Copynumber_meancol)$out)) 

#A boxplot can help us visualize the amount of outliers in our data. This is relevant as too many (extreme) outliers can
#have great impact on the results of our analysis and can change the outcome completely. They can easily affect the slope. 
#No outliers can be observed.

# (3) Density: Should be expect normality for drug sensitivity?

library(e1071)

par(mfrow=c(1, 2)) 

plot(density(lm_tab2$vorinostat), 
     main="Density Plot: Drug Sensitivity", 
     ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(lm_tab2$vorinostat), 2))
    )

polygon(density(lm_tab2$vorinostat), col="springgreen3")

#Skewness: 0.03 -> Plot is very slightly skewed to the right. 

plot(density(lm_tab2$BM_Copynumber_meancol), 
     main="Density Plot: Copynumber", 
     ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(lm_tab2$BM_Copynumber_meancol), 2))
    )

polygon(density(lm_tab2$BM_Copynumber_meancol), col="olivedrab1")

#Skewness: -0.29 -> Plot is slightly skewed to the left.


# (4) Correlation: what is the level level of linear dependence between the two variables?

cor(lm_tab2$vorinostat, lm_tab2$BM_Copynumber_meancol)

#[1] 0.06404028
#A good value for correlation lies close to 1 or -1, whilst the value 0 is undesirable. Values closer to 0 indicate
#that there is a weak relationship.  
#The value here is extremely low. 


#These analysis help us predcit whether a linear regression is or not the best model to describe our data. 
#Eventhough both the box plots (2) and the density plots (3) results could have been good indicators for a linear 
#relationship, the lack of a fitting straight line on the scatter plotc(1) and the low value in the result of the 
#correlation (4) indicate the opposite.

#Taking into consideration all results so far for this part, it is not unreasonable to predict that a linear
#regression will probably not be the best model to describe the relationships in our data.



###   4.2 Linear Regression                                                                                                 ####

## Linear Regression

reg2 <- lm(vorinostat ~ BM_Copynumber_meancol, data = lm_tab2)

## Details about the linear regression: what we need draw some conclusions

summary(reg2)

# Multiple R-squared:  0.004101 -> This indicates that only 0,4101% percent of the variation in the data (drug sensitivity)
#can be explained by the relationship between drug sensitivity and copynumber. In other words, there is a 0,4101% 
#variance reduction when we take the copynumber into account. 

# p-value: 0.6487 
#As the p-value for reg2 is significantly larger than 0.05 and R-squared tells us the copynumber only explains 0,4101% 
#of the variation in the data, it is safe to assume that there is no linear relationship between drug sensitivity and 
#copynumer, a.k.a copynumber cannot predict drug sensitivity. 

# More information about the fit (linear ecuation: y = y-intercept + slope * x) : 
confint(reg2)

###   4.3 Checking the normalization of residuals                                                                           ####

hist(reg2$residuals, 
     breaks = 20,
     xlab = "Residuals", 
     main = "Drug sensitivity vs copynumber: Histogram of the residuals")

# The data does NOT look normally distributed: As the residuals are not normally distributed, then the hypothesis that
#they are a random dataset, takes the value NO.
#This means that the linear regression model does not explain all trends in the dataset.

qqnorm(reg2$residuals)
qqline(reg2$residuals)


###   4.4 Visualization: Plots that describe the linear regression                                                          #### 


## Residual diagnostics: are the various assumptions that underpin linear regression reasonable for our data?

library(lattice)

xyplot(resid(reg2) ~ fitted(reg2),
       xlab = "Fitted Values",
       ylab = "Residuals",
       main = "Residual Diagnostic Plot",
       col = "hotpink",
       panel = function(x, y, ...)
       {
         panel.grid(h = -1, v = -1)
         panel.abline(h = 0)
         panel.xyplot(x, y, ...)
       }
      )


qqmath( ~ resid(reg2),
        xlab = "Theoretical Quantiles",
        ylab = "Residuals",
        abline(a=-2.8839, b=1))

#Here we expect to be able to draw a straight line that is fitting for most points. Although, this graphic seems promising,
#too many points would not be properly fit. 

# Visualization of regression 

par(mar = c(4, 4, 2, 2), mfrow = c(1, 2))
plot(reg2, which = c(1, 2))




####  5.  SIMPLE LINEAR REGRESSION WITH ALL GENES: Drug sensitivity with copynumber                                         ####

# Can we predict drug sensitivity using the copynumber data?
# How much of the variance of the data can be explained using the copynumber data?
## Data frames

CN_all = as.data.frame(CN_meancol)

DS = as.data.frame(drug_sensitivity)

## Table with drug sensitivity and doubling time per cell line

lm_tab_all2 = transform(merge(CN_all,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
lm_tab_all2 <- na.omit(lm_tab_all2)


###   5.1 Plots and visualization: Predicting how fit linear regression will be as a model to describe our data             #### 

## Checking for correlation

cor(lm_tab_all2)


## Ploting the data: can a linear relationship be observed? Should we have expected a high value for R-squared?

# Visualization of doubling time vs drug sensitivity


scatter.smooth(lm_tab_all2$vorinostat, 
               lm_tab_all2$CN_meancol, 
               col = "orangered",
               main = "Drug sensitivity & Copynumber Regression",
               xlab = "Drug sensitivity",
               ylab = "Copynumber",
               cex = 1.3,
               pch = 1)

#At first look, the points in the plot are so extremely scattered, that a linear relationship seems really unlikely. 
#The second problem that one can observe in this graphic, is that the line describing the behaviour
#is not straight. 


# (2) Box plot

par(mfrow=c(1, 2))
boxplot(lm_tab_all2$vorinostat, main="Drug sensitivity", sub=paste("Outlier rows: ", boxplot.stats(lm_tab_all2$vorinostat)$out)) 
boxplot(lm_tab_all2$CN_meancol, main="Copynumber", sub=paste("Outlier rows: ", boxplot.stats(lm_tab_all2$CN_meancol)$out)) 

#A boxplot can help us visualize the amount of outliers in our data. This is relevant as too many (extreme) outliers can
#have great impact on the results of our analysis and can change the outcome completely. They can easily affect the slope. 
#No outliers can be observed.

# (3) Density: Should be expect normality for drug sensitivity?

library(e1071)

par(mfrow=c(1, 2)) 

plot(density(lm_tab_all2$vorinostat), 
     main="Density Plot: Drug Sensitivity", 
     ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(lm_tab_all2$vorinostat), 2))
    )

polygon(density(lm_tab_all2$vorinostat), col="orangered")

#Skewness: 0.03 -> Plot is very slightly skewed to the right. 

plot(density(lm_tab_all2$CN_meancol), 
     main="Density Plot: Copynumber", 
     ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(lm_tab_all2$CN_meancol), 2))
    )

polygon(density(lm_tab_all2$CN_meancol), col="lightcoral")

#Skewness: -0.42 -> Plot is skewed to the left.


# (4) Correlation: what is the level level of linear dependence between the two variables?

cor(lm_tab2$vorinostat, lm_tab_all2$CN_meancol)

#[1] 0.1534729
#A good value for correlation lies close to 1 or -1, whilst the value 0 is undesirable. Values closer to 0 indicate
#that there is a weak relationship.  
#The value here is too low. 


#These analysis help us predcit whether a linear regression is or not the best model to describe our data. 
#Although no outliers were observed in the box plot results, all ther results seem to indicate that a linear 
#regression is most likely not the best model to describe our data. 


###   5.2 Linear Regression                                                                                                 ####

## Linear Regression

reg_all2 <- lm(vorinostat ~ CN_meancol, data = lm_tab_all2)

## Details about the linear regression: what we need draw some conclusions

summary(reg_all2)

# Multiple R-squared:  0.02355 -> This indicates that only 2.355% percent of the variation in the data (drug sensitivity)
#can be explained by the relationship between drug sensitivity and copynumber. In other words, there is a 2.355% 
#variance reduction when we take the copynumber into account. 

# p-value: 0.2726 
#As the p-value for reg2 is significantly larger than 0.05 and R-squared tells us the copynumber only explains 2.355% 
#of the variation in the data, it is safe to assume that there is no linear relationship between drug sensitivity and 
#copynumer, a.k.a copynumber cannot predict drug sensitivity. 

# More information about the fit (linear ecuation: y = y-intercept + slope * x) : 
confint(reg_all2)

###   5.3 Checking the normalization of residuals                                                                           ####


hist(reg_all2$residuals, 
     breaks = 20,
     xlab = "Residuals", 
     main = "Drug sensitivity vs copynumber: Histogram of the residuals")

# The data does NOT look normally distributed: As the residuals are not normally distributed, then the hypothesis that
#they are a random dataset, takes the value NO.
#This means that the linear regression model does not explain all trends in the dataset.

qqnorm(reg_all2$residuals)
qqline(reg_all2$residuals)

###   5.4 Visualization: Plots that describe the linear regression                                                          ####


# Residual diagnostics: are the various assumptions that underpin linear regression reasonable for our data?

library(lattice)

xyplot(resid(reg_all2) ~ fitted(reg_all2),
       xlab = "Fitted Values",
       ylab = "Residuals",
       main = "Residual Diagnostic Plot",
       col = "hotpink",
       panel = function(x, y, ...)
       {
         panel.grid(h = -1, v = -1)
         panel.abline(h = 0)
         panel.xyplot(x, y, ...)
       }
      )


qqmath( ~ resid(reg_all2),
        xlab = "Theoretical Quantiles",
        ylab = "Residuals",
        abline(a=-2.8839, b=1))

#Here we expect to be able to draw a straight line that is fitting for most points. Although, this graphic seems promising,
#too many points would not be properly fit. 

# Visualization of regression 

par(mar = c(4, 4, 2, 2), mfrow = c(1, 2))
plot(reg_all2, which = c(1, 2))



####  6.  MULTIPLE LINEAR REGRESSION WITH 100 BIOMARKERS: Drug sensitivity with doubling time and copynumber                ####


###   6.1 Plots and visualization: Predicting how fit linear regression will be as a model to describe our data             #### 
###   6.2 Linear Regression                                                                                                 ####

CN = as.data.frame(BM_Copynumber_meancol)

DT = as.data.frame(Doubling_Time)

DS = as.data.frame(drug_sensitivity)

lm_tab_m = transform(merge(CN, lm_tab,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_m <- na.omit(lm_tab_m)

reg_m <- lm(vorinostat ~ BM_Copynumber_meancol + Doubling_Time, data = lm_tab_m)

summary(reg_m)
###   6.3 Checking the normalization of residuals                                                                           ####
###   6.4 Visualization: Plots that describe the linear regression                                                          ####

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






##### PART 3 ##################################################################################################################
####  9.  Other models for regression

