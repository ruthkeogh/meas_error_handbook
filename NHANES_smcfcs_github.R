
#-----------------------------
#packages
#-----------------------------

#devtools::install_github("jwb133/smcfcs", ref="measerror")

library(survival)
library(flexsurv)
library(mitools)

#-----------------------------
#load data
#-----------------------------

mydata <- read.csv("nhanes.csv", header=TRUE)

#-----------------------------
#create complete-case data set
#-----------------------------

ccaData <- subset(mydata, (is.na(sbp1)==FALSE) & (is.na(smoke)==FALSE))
n <- dim(ccaData)[1]

head(ccaData)
# sbp1 sbp2 sex  age smoke diabetes d        t
# 5   0.5   NA   1 -0.7     1        0 0 12.00000
# 6  -0.1   NA   1 -1.0     1        0 0 21.16666
# 10  1.9   NA   1  2.0     0        0 0  3.25000
# 13  1.1  0.3   1  0.2     0        0 0 13.66666
# 15  0.8   NA   1  1.4     0        0 0  2.75000
# 16 -0.6   NA   1 -0.1     0        0 0  2.50000

#-----------------------------
#smcfcs using complete-case data - imputation using coxph
#-----------------------------

#create variable representing the true exposure: this is unobserved for all individuals
ccaData$sbp.true=NA

#predictor matrix: here just a matrix of zeros because no missing values
predmat=matrix(0,nrow=dim(ccaData)[2],ncol=dim(ccaData)[2])
colnames(predmat)=names(ccaData)
rownames(predmat)=names(ccaData)

#measurement error matrix
errmat=matrix(0,nrow=dim(ccaData)[2],ncol=dim(ccaData)[2])
colnames(errmat)=names(ccaData)
rownames(errmat)=names(ccaData)
errmat["sbp.true",c("sbp1","sbp2")]=1

#perform the multiple imputation
imp.1=smcfcs(originaldata=ccaData,smtype="coxph",smformula="Surv(t,d)~sbp.true+sex+age+smoke+diabetes",
           method=c("","","","","","","","","latnorm"),predictorMatrix=predmat,m=100,numit=100,
           rjlimit=1000,noisy=FALSE,errorProneMatrix=errmat)


#analysis using weibull model
impobj.1 <- imputationList(imp.1$impDatasets)
models.1<-with(impobj.1, flexsurvspline(Surv(t,d)~sbp.true+sex+age+smoke+diabetes,k=0, scale="hazard"))

#pool resulted across imputed data sets
summary(MIcombine(models.1))

#-----------------------------
#smcfcs using complete-case data - imputation using weibull
#-----------------------------

#create variable representing the true exposure: this is unobserved for all individuals
ccaData$sbp.true=NA

#predictor matrix: here just a matrix of zeros because no missing values
predmat=matrix(0,nrow=dim(ccaData)[2],ncol=dim(ccaData)[2])
colnames(predmat)=names(ccaData)
rownames(predmat)=names(ccaData)

#measurement error matrix
errmat=matrix(0,nrow=dim(ccaData)[2],ncol=dim(ccaData)[2])
colnames(errmat)=names(ccaData)
rownames(errmat)=names(ccaData)
errmat["sbp.true",c("sbp1","sbp2")]=1

#perform the multiple imputation
imp.1b=smcfcs(originaldata=ccaData,smtype="weibull",smformula="Surv(t,d)~sbp.true+sex+age+smoke+diabetes",
             method=c("","","","","","","","","latnorm"),predictorMatrix=predmat,m=100,numit=100,
             rjlimit=1000,noisy=FALSE,errorProneMatrix=errmat)


#analysis using weibull model
impobj.1b <- imputationList(imp.1b$impDatasets)
models.1b<-with(impobj.1b, flexsurvspline(Surv(t,d)~sbp.true+sex+age+smoke+diabetes,k=0, scale="hazard"))

#pool resulted across imputed data sets
summary(MIcombine(models.1b))

#-----------------------------
#smcfcs using all the data - imputation using coxph
#-----------------------------

#create variable representing the true exposure: this is unobserved for all individuals
mydata$sbp.true=NA

#predictor matrix: here we indicate which covariates are used in the imputation model for 'smoke'
predmat=matrix(0,nrow=dim(mydata)[2],ncol=dim(mydata)[2])
colnames(predmat)=names(mydata)
rownames(predmat)=names(mydata)
predmat["smoke",c("sbp.true","sex","age","diabetes")]=1

#measurement error matrix
errmat=matrix(0,nrow=dim(mydata)[2],ncol=dim(mydata)[2])
colnames(errmat)=names(mydata)
rownames(errmat)=names(mydata)
errmat["sbp.true",c("sbp1","sbp2")]=1

#perform the multiple imputation
imp.2=smcfcs(originaldata=mydata,smtype="coxph",smformula="Surv(t,d)~sbp.true+sex+age+smoke+diabetes",
           method=c("","","","","logreg","","","","latnorm"),predictorMatrix=predmat,m=100,numit=100,
           rjlimit=1000,noisy=FALSE,errorProneMatrix=errmat)

#analysis using weibull model
impobj.2 <- imputationList(imp.2$impDatasets)
models.2<-with(impobj.2, flexsurvspline(Surv(t,d)~sbp.true+sex+age+smoke+diabetes,k=0, scale="hazard"))

#pool resulted across imputed data sets
summary(MIcombine(models.2))

#-----------------------------
#smcfcs using all the data - imputation using weibull
#-----------------------------

#create variable representing the true exposure: this is unobserved for all individuals
mydata$sbp.true=NA

#predictor matrix: here we indicate which covariates are used in the imputation model for 'smoke'
predmat=matrix(0,nrow=dim(mydata)[2],ncol=dim(mydata)[2])
colnames(predmat)=names(mydata)
rownames(predmat)=names(mydata)
predmat["smoke",c("sbp.true","sex","age","diabetes")]=1

#measurement error matrix
errmat=matrix(0,nrow=dim(mydata)[2],ncol=dim(mydata)[2])
colnames(errmat)=names(mydata)
rownames(errmat)=names(mydata)
errmat["sbp.true",c("sbp1","sbp2")]=1

#perform the multiple imputation
imp.3=smcfcs(originaldata=mydata,smtype="weibull",smformula="Surv(t,d)~sbp.true+sex+age+smoke+diabetes",
             method=c("","","","","logreg","","","","latnorm"),predictorMatrix=predmat,m=100,numit=100,
             rjlimit=5000,noisy=FALSE,errorProneMatrix=errmat)

#analysis using weibull model
impobj.3 <- imputationList(imp.3$impDatasets)
models.3<-with(impobj.3, flexsurvspline(Surv(t,d)~sbp.true+sex+age+smoke+diabetes,k=0, scale="hazard"))

#pool resulted across imputed data sets
summary(MIcombine(models.3))

