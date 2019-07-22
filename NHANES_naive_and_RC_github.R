#-----------------------------
#packages
#-----------------------------

library(survival)
library(flexsurv)
library(lme4)
library(boot)

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
#naive weibull analysis 
#-----------------------------

naiveweibull=flexsurvspline(Surv(t,d)~sbp1+sex+age+smoke+diabetes,k=0, scale="hazard",data=ccaData)

#-----------------------------
#regression calibration analysis, ignoring missing data
#-----------------------------

n <-dim(ccaData)[1]

#regression calibration on complete cases
rc <- function(data, index) {
  locdata <- data[index,]
  locdata$id <- 1:n
  stackeddata <- rbind(locdata, locdata)
  stackeddata$sbp=NA
  stackeddata$sbp[1:n]<-locdata$sbp1
  stackeddata$sbp[(n+1):(2*n)]<-locdata$sbp2
  rc_mod <- lmer(sbp~1+sex+age+smoke+diabetes+(1|id), data=stackeddata)
  locdata$x_blup <- fitted(rc_mod)[1:n]
  rc_out_mod <- flexsurvspline(Surv(t,d)~x_blup+sex+age+smoke+diabetes,k=0, scale="hazard", data=locdata)
  return(rc_out_mod$coef)
}

rc(ccaData, 1:n)

set.seed(6723431)
boot <- boot(ccaData, rc, R=2000)
for (i in 1:7) {
  print(boot.ci(boot, index=i, type="perc"))
}
