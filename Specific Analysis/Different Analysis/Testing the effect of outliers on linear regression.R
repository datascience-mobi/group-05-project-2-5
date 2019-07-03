### Testing the effect of outliers on linear regression ###

# This document is a deeper exploration of some of the results obtained in the R Script "Linear Regressions",
#and while one does not need to read that document to understand this one or to be able to read this code,
#it is my advice to read that docuemnt first, as this one will focus on how outliers observed during the 
#"Linear Regressions" exploration affect our data and whether we should remove the outliers for better, 
#more accurate results or not.

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
###   2.1 Creating the necessary table (Cell lines with drug sensitivity)                                                   ####

## (1) Table 1: Selection of Doubling time from cellline_annotation

library(dplyr)

# Selecting the desired columns

Doubling_Time <- cellline_annotation %>%
  select(Cell_Line_Name, Doubling_Time)

# Changig the name of the "name" column to the names of the cell lines

row.names(Doubling_Time) <- Doubling_Time$Cell_Line_Name
Doubling_Time[1] <- NULL

Doubling_Time


####  3.  SIMPLE LINEAR REGRESSION WITH ALL GENES: Drug sensitivity with doubling time                                      ####

# Can we predict drug sensitivity using doubling time?
# How much of the variance of the data can be explained using the doubling time?

## Data frames
DT = as.data.frame(Doubling_Time)

DS = as.data.frame(drug_sensitivity)

## Table with drug sensitivity and doubling time per cell line

lm_tab = transform(merge(DT,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)


###   3.1 Boxplot: removing the outliers                                                                                    #### 

## Using the data from the R Script "Linear Regression", we know that outliers can be onserved in the boxplot
#for doubling time, therefore we focus on that data.

## Boxplot for doubling time
boxplot(lm_tab$Doubling_Time, main="Doubling time", sub=paste("Outlier rows: ", boxplot.stats(lm_tab$Doubling_Time)$out)) 

# What are the values of the outliers?

boxplot(lm_tab$Doubling_Time)$out
#[1] 66.8 79.5

outliers <- boxplot(mtcars$disp, plot=FALSE)$out

# Storing the values of the outliers in a vector

outliers_DT <- boxplot(lm_tab$Doubling_Time, plot=FALSE)$out

# Removing the outliers

lm_tab[which(lm_tab$Doubling_Time %in% outliers_DT),]
#We can see that the outliers correspond to those of the cell lines A498 (renal) and HOP-92 (lung).

lm_tab_out <- lm_tab[-which(lm_tab$Doubling_Time %in% outliers_DT),]

# checking our results

boxplot(lm_tab_out$Doubling_Time, main="Doubling time", sub=paste("Outlier rows: ", boxplot.stats(lm_tab_out$Doubling_Time)$out)) 

#There are no outliers now

###   3.2 Linear Regression                                                                                                 ####

# linear_regression = lm(predicted ~ predictor, data = dat) 

## Linear Regression without removing the outliers

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

## Linear Regression removing the outliers

reg_out <- lm(vorinostat ~ Doubling_Time, data = lm_tab_out)

## Details about the linear regression: what we need draw some conclusions

summary(reg_out)

# Multiple R-squared:  0.01767 -> This indicates that only 1,1767% percent of the variation in the data (drug sensitivity)
#can be explained by the relationship between drug sensitivity and doubling time. In other words, there is a 1.767% 
#variance reduction when we take the doubling time into account. 

# p-value: 0.3155 
#As the p-value for reg_out is significantly larger than 0.05 and R-squared tells us the doubling time only explains 1.767% 
#of the variation in the data, it is safe to assume that there is no linear relationship between drug sensitivity and 
#doubling time, a.k.a doubling time cannot predict drug sensitivity. 

##We can already observe this small change, the removal of outliers, has a significant impact on our results. 


# More information about the fit (linear ecuation: y = y-intercept + slope * x) : 

confint(reg1)


confint(reg_out)

#Small changes are observed.

###   3.3 Checking the normalization of residuals                                                                           ####

## Drug sensitivity vs doubling time: Histogram of the residuals

par(mfrow=c(1, 2))
hist(reg1$residuals, 
     breaks = 20, 
     xlab = "Residuals", 
     main = "Histogram with outliers")

hist(reg_out$residuals, 
     breaks = 20, 
     xlab = "Residuals", 
     main = "Histogram without outliers")

# The data in both histograms does not appear normalized and it is hard to say whether there
#is a case in which the data looks more normally distributed or not.
#This means that the linear regression models do not explain all trends in the dataset.

par(mfrow=c(1, 2))

qqnorm(reg1$residuals)
qqline(reg1$residuals, col = "red")

qqnorm(reg_out$residuals)
qqline(reg_out$residuals, col = "red")

# Small changes can be observed between values 1 and 2.

###   3.4 Visualization: Plots that describe the linear regression                                                          ####

## Residual diagnostics: are the various assumptions that underpin linear regression reasonable for our data?

library(lattice)

# 

par(mfrow=c(1, 2))

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



xyplot(resid(reg_out) ~ fitted(reg_out),
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

# The data for reg_out appear to be more spread.  

# 

par(mfrow=c(1, 2))

qqmath( ~ resid(reg1),
        xlab = "Theoretical Quantiles",
        ylab = "Residuals",
        abline(a=-0.01034, b=1))

qqmath( ~ resid(reg_out),
        xlab = "Theoretical Quantiles",
        ylab = "Residuals",
        abline(a=-0.01034, b=1))


# Visualization of regression 

par(mar = c(4, 4, 2, 2), mfrow = c(2, 2))
plot(reg1, which = c(1, 2))
plot(reg_out, which = c(1, 2))

# The line drawn accross the scatter plot changes when we remove the outliers.

# Comparing prediction and real values for drug sensitivity

par(mfrow=c(1, 2))

plot(lm_tab$vorinostat, reg1$fitted.values, pch = 20, col = "royalblue1", xlab = "Real values", ylab = "Predicted values")
abline(0, 1, col = "red")

plot(lm_tab_out$vorinostat, reg_out$fitted.values, pch = 20, col = "royalblue1", xlab = "Real values", ylab = "Predicted values")
abline(0, 1, col = "red")

# Changes in the red line can be observed_ change in the slope and intercept with the x-axis.



###   4.  General Conclusions                                                                                               ####