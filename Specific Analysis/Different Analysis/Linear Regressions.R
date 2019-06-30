

### Linear Regressions ###


## Do not run the entire r script at once. 
## Pending: Change log?
## No linear relationship: lowess, GAM
## Complete until 7

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

# linear_regression = lm(predicted ~ predictor, data = dat) 


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

# Comparing prediction and real values for drug sensitivity

plot(lm_tab$vorinostat, reg1$fitted.values, pch = 20, col = "royalblue1", xlab = "Real values", ylab = "Predicted values")
abline(0, 1, col = "red")


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

# Comparing prediction and real values for drug sensitivity

plot(lm_tab2$vorinostat, reg2$fitted.values, pch = 20, col = "springgreen3", xlab = "Real values", ylab = "Predicted values")
abline(0, 1, col = "red")



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
#As the p-value for reg_all2 is significantly larger than 0.05 and R-squared tells us the copynumber only explains 2.355% 
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

# Comparing prediction and real values for drug sensitivity

plot(lm_tab_all2$vorinostat, reg_all2$fitted.values, pch = 20, col = "lightcoral", xlab = "Real values", ylab = "Predicted values")
abline(0, 1, col = "royalblue1")



####  6.  MULTIPLE REGRESSION WITH 100 BIOMARKERS: Drug sensitivity with doubling time and copynumber                       ####

## Data frames

CN = as.data.frame(BM_Copynumber_meancol)

DT = as.data.frame(Doubling_Time)

DS = as.data.frame(drug_sensitivity)

## Table with drug sensitivity, copynumber and doubling time per cell line

lm_tab_m = transform(merge(CN, lm_tab,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_m <- na.omit(lm_tab_m)

names(lm_tab_m)[names(lm_tab_m) == "BM_Copynumber_meancol"] <- "Copynumber"
names(lm_tab_m)[names(lm_tab_m) == "vorinostat"] <- "Drug_Sensitivity"


###   6.1 Plots and visualization: Predicting how fit linear regression will be as a model to describe our data             #### 

## Checking for correlation

cor(lm_tab_m)

## Ploting the data: can a linear relationship be observed? Should we have expected a high value for R-squared?

# Visualization of doubling time vs drug sensitivity

# (1) Plot
plot(lm_tab_m , 
     pch=20 , 
     cex=1.5 , 
     col=rgb(0.5, 0.8, 0.9, 0.7)
    )

#There seems to be no linear relationship, because of how spread the points are.

# (2) Box plot

par(mfrow=c(1, 3))
boxplot(lm_tab_m$Drug_sensitivity, main="Drug sensitivity", sub=paste("Outlier rows: ", boxplot.stats(lm_tab_m$Drug_sensitivity)$out)) 
boxplot(lm_tab_m$Doubling_Time, main="Doubling time", sub=paste("Outlier rows: ", boxplot.stats(lm_tab_m$Doubling_Time)$out)) 
boxplot(lm_tab_m$Copynumber, main="Copynumber", sub=paste("Outlier rows: ", boxplot.stats(lm_tab_m$Copynumber)$out)) 

#A boxplot can help us visualize the amount of outliers in our data. This is relevant as too many (extreme) outliers can
#have great impact on the results of our analysis and can change the outcome completely. They can easily affect the slope. 
#Outliers can only be observed fo the boxplot of doubling time, which coincides with our previous results

# (3) Density: Should be expect normality for drug sensitivity?

library(e1071)

par(mfrow=c(1, 3)) 

plot(density(lm_tab_m$Drug_sensitivity), 
     main="Density Plot: Drug sensitivity", 
     ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(lm_tab_m$Drug_sensitivity), 2))
     )

polygon(density(lm_tab_m$Drug_sensitivity), col="darkorchid1")

#Skewness: 0.03 -> Plot is very slightly skewed to the right.

plot(density(lm_tab_m$Doubling_Time), 
     main="Density Plot: Doubling time", 
     ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(lm_tab_m$Doubling_Time), 2))
     )

polygon(density(lm_tab_m$Doubling_Time), col="skyblue1")

#Skewness: 1.02 -> Plot is slightly skewed to the left.

plot(density(lm_tab_m$Copynumber), 
     main="Density Plot: Copynumber", 
     ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(lm_tab_m$Copynumber), 2))
     )

polygon(density(lm_tab_m$Copynumber), col="olivedrab1")

#Skewness: -0.29 -> Plot is slightly skewed to the left.


#These analysis help us predcit whether a linear regression is or not the best model to describe our data. 
#Thanks to the outliers, considerably spread plots and skewed density plots, it is not unreasonable to predict that a 
#multiple regression with these parameterswill probably not be the best model to describe the relationships in our data.


###   6.2 Linear Regression                                                                                                 ####

## Linear Regression

reg_m <- lm(Drug_sensitivity ~ Copynumber + Doubling_Time, data = lm_tab_m)

## Details about the linear regression: what we need draw some conclusions

summary(reg_m)


# Multiple R-squared: 0.05949 -> This indicates that only 5.949% percent of the variation in the data (drug sensitivity)
#can be explained by the relationship between drug sensitivity, doubling time and copynumber. In other words, there is a 
#5.949% variance reduction when we take the both the doubling time and the copynumber into account. 


# p-value: 0.2158
#As the p-value for reg_m is significantly larger than 0.05 and R-squared tells us the copynumber only explains 2.355% 
#of the variation in the data, it is safe to assume that there is no linear relationship between drug sensitivity and 
#copynumer, a.k.a copynumber cannot predict drug sensitivity. 


# Coefficients:
#                   Estimate Std.    Error     t-value   Pr(>|t|)    
#   (Intercept)      6.139792       0.115579    53.122    <2e-16 ***
#   Copynumber       0.373463       0.834851    0.447     0.6566    
#   Doubling_Time   -0.005105       0.002975   -1.716     0.0924

#F-statistic (multiple regression) : 1.581
#F-statistic (drug sensitivity vs doubling time): 2.609 
#F-statistic (drug sensitivity vs copynumber):  0.21

#The coefficients show that better results are yielded when using doubling time alone to predict drug sensitivity,
#than when using both doubling time and copynumber.


## RMSE: A measure of how accurately the model predicts the response
n = nrow(lm_tab_m)
rmse = sqrt(1/n * sum(reg_m$residuals^2))
rmse

#[1] 0.2877142
#Low values indicate a better fit. This result is surprising as the R-squared and p-value tell us that these variables
#do not contribute significantly to the variance of the data. 

# More information about the fit (linear ecuation: y = y-intercept + slope * x) : 
confint(reg_m)

###   6.3 Checking the normalization of residuals                                                                           ####

hist(reg_m$residuals, 
     breaks = 20,
     xlab = "Residuals", 
     main = "Drug sensitivity vs copynumber and doubling time: Residuals histogram")

# The data does NOT look normally distributed: As the residuals are not normally distributed, then the hypothesis that
#they are a random dataset, takes the value NO.
#This means that the linear regression model does not explain all trends in the dataset.

qqnorm(reg_m$residuals)
qqline(reg_m$residuals)

###   6.4 Visualization: Plots that describe the linear regression                                                          ####

## Residual diagnostics: are the various assumptions that underpin linear regression reasonable for our data?

library(lattice)

xyplot(resid(reg_m) ~ fitted(reg_m),
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


qqmath( ~ resid(reg_m),
        xlab = "Theoretical Quantiles",
        ylab = "Residuals",
        abline(a=-2.8839, b=1))

#Here we expect to be able to draw a straight line that is fitting for most points. Although, this graphic seems promising,
#too many points would not be properly fit. 

# Visualization of regression 

par(mar = c(4, 4, 2, 2), mfrow = c(1, 2))
plot(reg_m, which = c(1, 2))


# Comparing prediction and real values for drug sensitivity

plot(lm_tab_m$Drug_sensitivity, reg_m$fitted.values, pch = 20, col = "orchid1", xlab = "Real values", ylab = "Predicted values")
abline(0, 1, col = "orange1")

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

