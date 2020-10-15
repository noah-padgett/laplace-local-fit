
library(coda)
library(mvtnorm)
library(lavaan)
library(data.table)
library(dplyr)
library(tcltk)

wd <- getwd()
source(paste0(wd, "/code/utility_functions.R"))


# specify population model
population.model <- ' 
  f1 =~ 1*y1 + 0.8*y2 + 0.8*y3 + 1.2*y4 + 0.02*y5 + -0.05*y6 + 0.05*y7 + 0.01*y8 + -0.02*y9 + 0.01*y10 + 0.03*y11 + 0.01*y12
  f2 =~ 1*y5 + 1.1*y6 + 0.8*y7 + 0.9*y8 + 0.01*y1 + -0.2*y2 + 0.02*y3 + 0.02*y4 + -0.02*y9 + 0.01*y10 + 0.03*y11 + 0.01*y12
  f3 =~ 1*y9 + 0.8*y10 + 1.3*y11 + 0.8*y12 + 0.01*y1 + -0.2*y2 + 0.02*y3 + 0.02*y4 + 0.02*y5 + -0.05*y6 + 0.05*y7 + 0.01*y8
  
  # Factor (co)variances
  f1 ~~ 0.5*f2 + 0.4*f3
  f2 ~~ 0.3*f3
  
  # residual covariances
  y1 ~~ 0.01*y2
  y2 ~~ 0.3*y3
  y4 ~~ 0.3*y7
'


nItems <- 12
nFactors <- 3
loadVal <- matrix(nrow=nItems, ncol=nFactors)
loadVal[1:nItems, 1:nFactors] <- runif(nItems*nFactors, -0.05, 0.05)
loadVal[1:4,  1] <- runif(4, 0.6, 0.8)
loadVal[5:8,  2] <- runif(4, 0.6, 0.8)
loadVal[9:12, 3] <- runif(4, 0.6, 0.8)
load <- matrix(nrow=nItems, ncol=nFactors)
LY <- bind(load, loadVal)

# latent variable covariance matrix
lvc <- matrix(nrow=nFactors, ncol=nFactors)
diag(lvc) <- 1
lvcVal <- matrix(nrow=nFactors, ncol=nFactors)
diag(lvcVal) <- 1
lvcVal[is.na(lvcVal)] <- 0.3
RPS <- binds(lvc, lvcVal)

# error covariance matrix
ecov <- matrix(0,nrow=nItems, ncol=nItems)
diag(ecov) <- NA
ecov[1,4] <- ecov[4,1] <- NA
ecov[5,6] <- ecov[6,5] <- NA
ecovVal <- matrix(nrow=nItems, ncol=nItems)
diag(ecovVal) <- 1
ecovVal[1,4] <- ecovVal[4,1] <- 0.4
ecovVal[5,6] <- ecovVal[6,5] <- 0.4

RTE <- binds(ecov, ecovVal)

CFA.MODEL <- model(LY=LAMBDA, RPS=PHI, RTE=PSI, modelType = "CFA")


# generate data
set.seed(1234)
myData <- simulateData(population.model, sample.nobs=300L)

# population moments
fitted(sem(population.model))

# sample moments
round(cov(myData), 3)
round(colMeans(myData), 3)

# fit model
myModel <- ' 
  f1 =~ y1 + y2 + y3 + y4
  f2 =~ y5 + y6 + y7 + y8 
  f3 =~ y9 + y10 + y11 + y12 
  
  # Factor covariances
  f1 ~~ f2 + f3
  f2 ~~ f3
  
  # residual covariances
  y2 ~~ y3
  y4 ~~ y7
  
'
fit <- cfa(myModel, data=myData)
summary(fit)



# Monte Carlo simulation using simsem
library(simsem)


Output <- sim(100, model=myModel, n=200, generate=population.model, lavaanfun = "cfa")

mean(Output@stdCoef$`y2~~y3` >= 0.25)
summary(Output@stdCoef$`y2~~y3`)

mean(Output@stdCoef$`y4~~y7` >= 0.25)
summary(Output@stdCoef$`y4~~y7`)
