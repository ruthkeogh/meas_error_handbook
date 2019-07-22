
#-----------------------------
#packages
#-----------------------------

library(survival)
library(flexsurv)
library(lme4)
library(boot)
library(MASS)
library(rjags)

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

naiveweibull <- survreg(Surv(t,d)~sbp1+sex+age+smoke+diabetes,dist="weibull", data=ccaData)

#-----------------------------
#Bayes analysis, allowing for measurement error but ignoring missing data
#-----------------------------

set.seed(87893)
nChains <- 5

#construct vectors to model as censored weibull
t <- ccaData$t
t[ccaData$d==0] <- NA
c <- ccaData$t
iscensored <- 1-ccaData$d

#fit measurement error model to get some overdispersed initial values
wideccaData <- ccaData
wideccaData$id <- 1:n
wideccaData$meas <- 1
longData <- rbind(wideccaData, wideccaData)
longData$sbp1[(n+1):(2*n)] <- ccaData$sbp2
longData$meas[(n+1):(2*n)] <- 2
rc_mod <- lmer(sbp1~1+sex+age+smoke+diabetes+factor(meas)+(1|id), data=longData)

#initial values

#JAGS treats the event times of those subjects who were censored as missing data
#we therefore need to specify sensible initial values for these
jags.inits <- vector("list",nChains)
tInits <- array(0, dim=c(nChains,n))
logscaleinits <- array(0, dim=c(nChains, 6))
rinits <- array(0, dim=c(nChains))

for (i in 1:nChains) {
  tInits[i,] <- rep(NA,n)
  tInits[i,iscensored==1] <- ccaData$t[iscensored==1]+0.5*runif(sum(iscensored))
  logscaleinits[i,] <- mvrnorm(1, coef(naiveweibull), 2*vcov(naiveweibull)[1:6,1:6])
  newlogscale <- rnorm(1, log(naiveweibull$scale), 2*vcov(naiveweibull)[7,7]^0.5)
  rinits[i] <- 1/exp(newlogscale)
  jags.inits[[i]] <- list(t=tInits[i,], log.scale=logscaleinits[i,], r=rinits[i])
}

gammaInits <- array(0, dim=c(nChains,5))
tauxInits <- array(0, dim=nChains)
tauuInits <- array(0, dim=nChains)
for (i in 1:nChains) {
  draw <- mvrnorm(1, fixef(rc_mod), 2*vcov(rc_mod))
  gammaInits[i,] <- draw[1:5]
  tauxInits[i] <- ((6-i)/2.5)/(VarCorr(rc_mod)$id[1]^2)
  tauuInits[i] <- (i/2.5)/(sigma(rc_mod)^2)
  jags.inits[[i]] <- list(t=tInits[i,], log.scale=logscaleinits[i,], r=rinits[i], gamma=gammaInits[i,], taux=tauxInits[i], tauu=tauuInits[i])
}


gammamean <- rep(0,5)
gammaprec <- 0.0001*diag(5)

jags.data<-list(n = n, t=t, c=c, is.censored=iscensored,
                sbp1 = ccaData$sbp1, sbp2 = ccaData$sbp2, sex=ccaData$sex, age=ccaData$age, smoke=ccaData$smoke, diabetes=ccaData$diabetes,
                "gammamean"=gammamean, "gammaprec"=gammaprec,"taux_alpha"=0.5, "taux_beta"=0.5,"tauu_alpha"=0.5, "tauu_beta"=0.5)

jagsmodel2 <- jags.model(data=jags.data, file="weibull_adj_cca.bug",n.chains=nChains, inits=jags.inits)
burnin2 <- coda.samples(jagsmodel2, variable.names=c("beta","r","sigmau","sigmax"), n.iter=5000, thin=1)
gelman.diag(burnin2)
summary(burnin2)
mainsample2 <- coda.samples(jagsmodel2, variable.names=c("beta","r","sigmau","sigmax"), n.iter=5000, thin=1)
gelman.diag(mainsample2)
summary(mainsample2)



