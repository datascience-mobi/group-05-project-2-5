################## Creat Vorinostat data ####################

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

################################################################################################

################################################################################################

##### check normalization ####
# see alson in Specific Analysis T test

#Anyway: Check normality of the FC-values via qqplot:
FC <- TreatedVorinostat - UntreatedVorinostat
qqnorm(FC)
qqline(FC, col= "red")
#We see, that we have a heavy tailed distribution. Nevertheless, we assume normality.


#For normalization, a package needs to be installed: install.packages("BBmisc")
library(BBmisc)
FCnorm <- normalize(FC, method= "scale")
#Check, if data look more normal distributed now:
qqnorm(FCnorm)
qqline(FCnorm, col= "red")
#-> more linear than before

# same check for treated and untreated 
qqnorm(TreatedVorinostat)
qqline(TreatedVorinostat, col= "red")

qqnorm(UntreatedVorinostat)
qqline(UntreatedVorinostat, col= "red")