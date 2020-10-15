# ============================ #
# Bayesian CFA Approximation
#   with Laplace method
#   Generalized version
# ============================ #
# Created by: R. Noah Padgett
# Created on: 2020-02-18
# 
# Laste Editted: 2020-02-18
# ============================ #


library(coda)
library(mvtnorm)
library(lavaan)
library(data.table)
library(dplyr)
library(tcltk)

wd <- getwd()
source(paste0(wd, "/code/utility_functions.R"))

# ========================================== #
# ========================================== #
#   function: laplace_approx()
# ========================================== #
# use: runs the laplace approximate for a 
#       given parameter
#
# arguments:
# model      - list of model components (lambda...)
# inits      - initial values
# no.samples - 
#
laplace_approx <- function(model, inits, no.samples, scale.cov, ...) {
    #fit <- nlminb(inits, model,
  #               control=list(eval.max=20000, iter.max=10000,
  #                            abs.tol=2.2e-15, rel.tol=1e-10,
  #                            x.tol=1.5e-8,xf.tol=2.2e-14), ...)
  
  fit <- optim(inits, model, control = list(fnscale = -1), hessian = TRUE, ...)
  param_mean <- fit$par # numerical deriv
  # compute hess at param_mean
  #hess <- numDeriv::hessian(model, param_mean, ...)
  #param_cov_mat <- solve(-hess)
  param_cov_mat <- solve(-fit$hessian)
  
  # scaled covariance matrix (artifically inflate uncertainty)
  A <- diag(scale.cov, nrow=nrow(param_cov_mat), ncol=ncol(param_cov_mat))
  param_cov_mat <- A%*%param_cov_mat%*%t(A)

  mcmc(rmvnorm(no.samples, param_mean, param_cov_mat))
}


# ========================================== #
# ========================================== #
#   function: get_prior_dens()
# ========================================== #
# use: gets the appropriate prior for the 
#       parameter of interest
#
get_prior_dens <- function(pvalue, pname,...){
  if(pname %like% 'lambda'){
    out <- dnorm(pvalue, 0, 1, log=T)
  }
  if(pname %like% 'dphi'){
    out <- dgamma(pvalue, 1, 0.5, log=T)
  }
  if(pname %like% 'odphi'){
    out <- dnorm(pvalue, 0, 1, log=T)
  }
  if(pname %like% 'dpsi'){
    out <- dgamma(pvalue, 1, 0.5, log=T)
  }
  if(pname %like% 'odpsi'){
    out <- dnorm(pvalue, 0, 1, log=T)
  }
  if(pname %like% 'eta'){
    out <- dnorm(pvalue, 0, 10, log=T)
  }
  if(pname %like% 'tau'){
    out <- dnorm(pvalue, 0, 32, log=T)
  }
  return(out)
}

# ========================================== #
# ========================================== #
#   function: get_log_post()
# ========================================== #
# use: uses the model, parameters, and data to
#       to calculate log posterior
#
# arguments:
# p        - names vector of parameters
# X        - data frame of raw data
# cfaModel - list of model components
#
get_log_post <- function(p, X, cfaModel,...) {
  out <- use_cfa_model(p, X, cfaModel)
  log_lik <- sum(apply(X, 1, dmvnorm, mean=out[['tau']], sigma=out[['Sigma']], log=T))
  
  log_prior<-0
  if(length(p)==1){
    log_prior <- get_prior_dens(p, names(p))
  } else {
    i <- 1
    for(i in 1:length(p)){
      log_prior <- log_prior + get_prior_dens(p[i], names(p)[i])
    }
  }
  log_post <- log_lik + log_prior
  log_post
}


# ========================================== #
# ========================================== #
#   function: laplace_local_fit()
# ========================================== #
# use: uses the fittes lavaan object to run
#       the proposed method
#
# arguments:
# fit       - fitted lavaan model
# cut.load  - cutoff for value of loading to care about default = 0.3 
# cut.cov   - cutoff for value of covariances to care about default = 0.1
# opt       - list of parameters to pass to interior functions
# sum.print - logical indicator of whether to print the summary table upon completion
# counter   - logical indicator of whether to print out a (.) after each
#               parameter is completed
#
laplace_local_fit <- function(fit, cut.load = 0.3, cut.cov = 0.1,
                              opt=list(scale.cov=1, no.samples=1000),
                              sum.print=F, pb=T,...){
  
  extractedLavaan <- lavMatrixRepresentation(partable(fit))
  
  factNames <- unique(extractedLavaan[extractedLavaan[,"mat"]=="lambda", "lhs"])
  varNames <- unique(extractedLavaan[extractedLavaan[,"mat"]=="lambda", "rhs"])
  # extract factor loading matrix
  lambda <- extractedLavaan[ extractedLavaan$mat == "lambda" ,]
  lambda <- convert2matrix(lambda$row, lambda$col, lambda$est)
  colnames(lambda) <- factNames
  rownames(lambda) <- varNames
  # extract factor covariance matrix
  phi <- extractedLavaan[ extractedLavaan$mat == "psi" ,]
  phi <- convert2matrix(phi[,'row'], phi[,'col'], phi[,'est'])
  phi <- up2full(phi)
  colnames(phi) <- rownames(phi) <- factNames
  # extract error covariance matrix
  psi <- extractedLavaan[ extractedLavaan$mat == "theta" ,]
  psi <- convert2matrix(psi[,'row'], psi[,'col'], psi[,'est'])
  psi[upper.tri(psi)] <- 0
  colnames(psi) <- rownames(psi) <- varNames
  
  
  # need to create list of all NA parameters in the above matrices
  lamList <- as.matrix(which(is.na(lambda), arr.ind = T))
  il <- nrow(lamList)
  phiList <- as.matrix(which(is.na(phi), arr.ind = T))
  ip <- il + nrow(phiList)
  psiList <- as.matrix(which(is.na(psi), arr.ind = T))
  it <- ip + nrow(psiList)
  modList <- rbind(lamList, phiList, psiList)
  # number of variables
  # create names for each condition
  vnlamList <- lamList
  vnlamList[,2] <- paste0(factor(vnlamList[,2], levels = order(unique(vnlamList[,2])),labels=factNames))
  vnlamList[,1] <- rownames(lamList)
  vnlamList[,2] <- paste0(vnlamList[,2],"=~",vnlamList[,1])
  vnphiList <- phiList
  if(nrow(phiList)>0){
    vnphiList[,1] <- paste0(factor(phiList[,1], levels = order(unique(vnphiList[,1])),labels=factNames))
    vnphiList[,2] <- paste0(factor(phiList[,2], levels = order(unique(phiList[,2])),labels=factNames))
  }
  vnpsiList <- psiList
  vnpsiList[,1] <- rownames(psiList)
  vnpsiList[,2] <- paste0(vnpsiList[,1],"~~x", psiList[,2])
  nameList <- rbind(vnlamList, vnphiList, vnpsiList)
  
  # set up function that runs the approximation
  run_approximation <- function(x, lambdaMod = lambda, phiMod = phi, psiMod = psi, dots,...){
    # from here down, we need to loop around it.
    # x <- modList[i,]
    # lambdaMod <- lambda
    # phiMod <- phi
    # psiMod <- psi
    #data <- data
    # do we need to update lamda?
    if(i <= il){
      Q <- lambda
      Q[is.na(Q)] <- 0
      Q[x[1], x[2]] <- NA
      lambdaMod <- Q
    } else {
      Q <- lambda
      Q[is.na(Q)] <- 0
      lambdaMod <- Q
    }
    
    # update phi?
    if(i > il & i <= ip){
      Q <- phi
      Q[is.na(Q)] <- 0
      Q[x[1], x[2]] <- NA
      phiMod <- Q
    } else {
      Q <- phi
      Q[is.na(Q)] <- 0
      phiMod <- Q
    }
    
    # update psi?
    if(i > ip){
      Q <- psi
      Q[is.na(Q)] <- 0
      Q[x[1], x[2]] <- NA
      psiMod <- Q
    } else {
      Q <- psi
      Q[is.na(Q)] <- 0
      psiMod <- Q
    }
    
    # combine into a single list
    cfaModel <- list(lambdaMod, phiMod, psiMod) #, tauMod, etaMod
    
    #print(cfaModel)
    # get starting values
    inits <- get_starting_values(cfaModel) 
    
    # run the approximatio
    samples <- laplace_approx(get_log_post, inits, X=data, scale.cov = dots[[1]],
                              no.samples=dots[[2]], cfaModel=cfaModel)
    return(samples)
  }
  
  # iterate around this function
  i <- 1
  fitResults <- matrix(nrow=opt[[2]], ncol=it)
    # progress bar
    progress_bar <- txtProgressBar(min = 0, max = it, style = 3)
  for(i in 1:it){
    fitResults[,i] <- run_approximation(modList[i,], dots=opt)
    if(pb == T) setTxtProgressBar(progress_bar, i)
  }
  colnames(fitResults) <- nameList[,2, drop=T]
  # now, compute and format summary statistics
  sumResults <- data.frame(matrix(nrow=ncol(fitResults), ncol=9))
  colnames(sumResults) <- c("Parameter","Prob", "mean", "sd", "p0.025", "p0.25", "p0.5", "p0.75", "p0.975")
  sumResults[,1] <- colnames(fitResults)
  
  sumResults[,3:9] <- t(apply(fitResults, 2, function(x){
    c(mean(x, na.rm=T), sd(x, na.rm=T),
      quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=T))
  }))
  
  # compute probability of meaningfulness
  # depends on parameter
  # cut.load = 0.3
  # cut.cov = 0.1
  p <- colnames(fitResults)
  for(i in 1:ncol(fitResults)){
    x <- fitResults[,i, drop=T]
    if(p[i] %like% "=~"){
      pv <- mean(ifelse(abs(x) > 0.3, 1, 0))
    }
    if(p[i] %like% "~~"){
      pv <- mean(ifelse(abs(x) > 0.1, 1, 0))
    }
    sumResults[i, 2] <- pv
  }
  sumResults <- arrange(sumResults, desc(Prob))

  sumResults[,2:9] <- round(sumResults[,2:9], 3)
  cat("\n")
  if(sum.print==T) print(sumResults, row.names = FALSE)
  
  # convertto data.frame
  fitResults <- as.data.frame(fitResults)
  out <- list(fitResults, sumResults)
  names(out) <- c("All Results", "Summary")
  
  return(out)
}



HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939)
summary(fit, fit.measures = TRUE)

data <- lavaan::HolzingerSwineford1939[, paste0('x',1:9)]


lfit <- laplace_local_fit(fit, data=data, cut.load = 0.4, cut.cov = 0.2,
                          opt=list(scale.cov=21, no.samples=100))
print(lfit$Summary)

plot(density(lfit$`All Results`$`textual=~x2`))
plot(density(lfit$`All Results`$`x5~~x2`))

