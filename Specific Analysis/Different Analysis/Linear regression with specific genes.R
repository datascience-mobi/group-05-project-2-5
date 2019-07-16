
### Linear Regression with top 3 Biomarkers in table BM_copynumber
###############

DHRS2 <- BM_copynumber[-c(2:90),] 


### Data frames

DHRS2 = as.data.frame(t(DHRS2))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_DHRS2 = transform(merge(DHRS2,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_DHRS2 <- na.omit(lm_tab_DHRS2)


reg_DHRS2 <- lm(vorinostat ~ DHRS2, data = lm_tab_DHRS2)


### Details about the linear regression: what we need draw some conclusions

summary(reg_DHRS2)


##############



CDKN1A <- BM_copynumber[-c(3:90),] 
CDKN1A <- CDKN1A[-c(1),] 



### Data frames

CDKN1A = as.data.frame(t(CDKN1A))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_CDKN1A = transform(merge(CDKN1A,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_CDKN1A <- na.omit(lm_tab_CDKN1A)


reg_CDKN1A <- lm(vorinostat ~ CDKN1A, data = lm_tab_CDKN1A)


### Details about the linear regression: what we need draw some conclusions

summary(reg_CDKN1A)


#################


CITED2 <- BM_copynumber[-c(4:90),] 
CITED2 <- CITED2[-c(1:2),] 



### Data frames

CITED2 = as.data.frame(t(CITED2))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_CITED2 = transform(merge(CITED2,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_CITED2 <- na.omit(lm_tab_CITED2)


reg_CITED2 <- lm(vorinostat ~ CITED2, data = lm_tab_CITED2)


### Details about the linear regression: what we need draw some conclusions

summary(reg_CITED2)



### Linear Regression with the 3 non-biomarker genes
###############

#################

CLPS <- Copynumber["CLPS",]


### Data frames

CLPS = as.data.frame(t(CLPS))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_CLPS = transform(merge(CLPS,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_CLPS <- na.omit(lm_tab_CLPS)


reg_CLPS <- lm(vorinostat ~ CLPS, data = lm_tab_CLPS)


### Details about the linear regression: what we need draw some conclusions

summary(reg_CLPS)


#################
#################

UPK2 <- Copynumber["UPK2",]


### Data frames

UPK2 = as.data.frame(t(UPK2))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_UPK2 = transform(merge(UPK2,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_UPK2 <- na.omit(lm_tab_UPK2)


reg_UPK2 <- lm(vorinostat ~ UPK2, data = lm_tab_UPK2)


### Details about the linear regression: what we need draw some conclusions

summary(reg_UPK2)


#################

#Largest Copynumber values


#We work with mean of the rows because we only want to compare the genes 
Copynumber_av= abs(Copynumber)
Copynumber_mean= rowMeans(Copynumber_av)

Copynumber_sav <- sort(Copynumber_mean, decreasing = TRUE)
Copynumber_sav <- as.matrix(Copynumer_sav)

#Top 10
CN_top10 = Copynumber_sav[1:10,]
CN_top10

#Last 10
CN_last10 = Copynumber_sav[23306:23316,]
CN_last10


### Compare linaer regression values for top 3 and last 3


###TOP 1

DAZ2 <- Copynumber["DAZ2",]




### Data frames

DAZ2 = as.data.frame(t(DAZ2))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_DAZ2 = transform(merge(DAZ2,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_DAZ2 <- na.omit(lm_tab_DAZ2)


reg_DAZ2 <- lm(vorinostat ~ DAZ2, data = lm_tab_DAZ2)


### Details about the linear regression: what we need draw some conclusions

summary(reg_DAZ2)


###TOP 2

UTY <- Copynumber["UTY",]

### Data frames

UTY = as.data.frame(t(UTY))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_UTY = transform(merge(UTY,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_UTY <- na.omit(lm_tab_UTY)


reg_UTY <- lm(vorinostat ~ UTY, data = lm_tab_UTY)


### Details about the linear regression: what we need draw some conclusions

summary(reg_UTY)

###TOP 3

DAZ1 <- Copynumber["DAZ1",]

### Data frames

DAZ1 = as.data.frame(t(DAZ1))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_DAZ1 = transform(merge(DAZ1,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_DAZ1 <- na.omit(lm_tab_DAZ1)


reg_DAZ1 <- lm(vorinostat ~ DAZ1, data = lm_tab_DAZ1)


### Details about the linear regression: what we need draw some conclusions

summary(reg_DAZ1)


###LAST 3

RBM18 <- Copynumber["RBM18",]

### Data frames

RBM18 = as.data.frame(t(RBM18))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_RBM18 = transform(merge(RBM18,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_RBM18 <- na.omit(lm_tab_RBM18)


reg_RBM18 <- lm(vorinostat ~ RBM18, data = lm_tab_RBM18)


### Details about the linear regression: what we need draw some conclusions

summary(reg_RBM18)


###LAST 2

MIR3974 <- Copynumber["MIR3974",]

### Data frames

MIR3974 = as.data.frame(t(MIR3974))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_MIR3974 = transform(merge(MIR3974,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_MIR3974 <- na.omit(lm_tab_MIR3974)


reg_MIR3974 <- lm(vorinostat ~ MIR3974, data = lm_tab_MIR3974)


### Details about the linear regression: what we need draw some conclusions

summary(reg_MIR3974)



###LAST 1

SKP1P2 <- Copynumber["SKP1P2",]

### Data frames

SKP1P2 = as.data.frame(t(SKP1P2))

DS = as.data.frame(drug_sensitivity)


### Table with drug sensitivity and doubling time per cell line

lm_tab_SKP1P2 = transform(merge(SKP1P2,DS,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

lm_tab_SKP1P2 <- na.omit(lm_tab_SKP1P2)


reg_SKP1P2 <- lm(vorinostat ~ SKP1P2, data = lm_tab_SKP1P2)


### Details about the linear regression: what we need draw some conclusions

summary(reg_SKP1P2)






