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
# sample.data - data frame of raw data
# cfa.model - list of model components
#
get_log_post <- function(p, sample.data, cfa.model,...) {
  
  out <- use_cfa_model(p, cov(sample.data), cfa.model)
  
  log_lik <- sum(apply(sample.data, 1, dmvnorm,
                       mean=out[['tau']],
                       sigma=out[['Sigma']], log=T))
  
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
#   function: use_cfa_model()
# ========================================== #
# use: take in parameters, data, and model to 
#         obtain the log-likelihood
#
# arguments:
# theta - vector of parameters being optimized
# sample.cov - samplecovariance matrix
# cfa.model - list of model parameters
use_cfa_model <- function(theta, sample.cov, cfa.model,...){
  # Compue sample statistics
  p<-ncol(sample.cov)
  S<-sample.cov
  
  # unpack model
  lambda <- cfa.model[[1]]
  phi <- cfa.model[[2]]
  psi <- cfa.model[[3]]
  #tau <- cfaModel[[4]]
  #eta <- cfaModel[[5]]
  
  # number factor loadings
  lam.num <- length(which(is.na(lambda)))
  lambda[which(is.na(lambda))] <- theta[1:lam.num]
  nF = ncol(lambda)
  # number elements in factor (co)variance matrix
  phi.num <- length(which(is.na(phi)))
  dphi.num <- sum(is.na(diag(phi))==T)
  odphi.num <- sum(is.na(phi[lower.tri(phi)])==T)
  if(phi.num > 0){
    if(dphi.num == 0){
      phi[which(is.na(phi))] <- theta[(lam.num+1):(lam.num+phi.num)]
    } else {
      diag(phi) <- theta[(lam.num+1):(lam.num+dphi.num)]
      phi[which(is.na(phi))] <- theta[(lam.num+dphi.num+1):(lam.num+phi.num)]
    }
  }
  phi <- low2full(phi) # map lower to upper
  
  # number elements in error (co)variance matrix
  psi.num <- length(which(is.na(psi)))
  dpsi.num <- sum(is.na(diag(psi))==T)
  odpsi.num <- sum(is.na(psi[lower.tri(psi)])==T)
  if(psi.num > 0){
    if(dpsi.num == 0){
      psi[which(is.na(psi))] <- theta[(lam.num+1):(lam.num+psi.num)]
    } else {
      diag(psi) <- theta[(lam.num+1):(lam.num+dpsi.num)]
      psi[which(is.na(psi))] <- theta[(lam.num+dpsi.num+1):(lam.num+psi.num)]
    }
  }
  psi <- low2full(psi)
  # number of factor scores
  #eta.num <- length(eta)
  #eta <- matrix(theta[(lam.num+phi.num+psi.num+tau.num+1):(lam.num+phi.num+psi.num+tau.num+eta.num)],
  #              nrow=nF)
  # mean center eta
  #for(i in 1:nF){
  #  eta[i, ] <- eta[i,] - mean(eta[,i])
  #}
  
  # # number of intercepts
  # tau.num <- length(tau)
  # tau <- matrix(theta[(lam.num+phi.num+psi.num+1):(lam.num+phi.num+psi.num+tau.num)], ncol=1)
  # tau <- repeat_col(tau, ncol(eta))
  
  # compute model observed outcomes
  #Y <- tau + lambda%*%eta
  tau <- numeric(p)
  # compute model implied (co)variance matrix
  Sigma<-lambda%*%phi%*%(t(lambda)) + psi
  
  #return fit value 
  out <- list(Sigma, lambda, phi, psi, tau)
  names(out) <- c('Sigma', 'lambda', 'phi', 'psi', 'tau')
  return(out)
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
# standardized - logical for whether to standardize
# cut.load  - cutoff for value of loading to care about default = 0.3 
# cut.cov   - cutoff for value of covariances to care about default = 0.1
# opt       - list of parameters to pass to interior functions
# sum.print - logical indicator of whether to print the summary table upon completion
# counter   - logical indicator of whether to print out a (.) after each
#               parameter is completed
#
laplace_local_fit <- function(fit, cut.load = 0.3, cut.cov = 0.1, standardize=T,
                              opt=list(scale.cov=1, no.samples=1000),
                              all.parameters=F,
                              sum.print=F, pb=T,...){
  
  # Observed Data
  sampleData <- fit@Data@X[[1]]
  # sample covariance matrix
  sampleCov <- fit@SampleStats@cov[[1]]
  
  # extract model
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
  
  if(all.parameters == T){
    lambdaA <- lambda
    phiA <- phi
    psiA <- psi
    
    lambdaA[!is.na(lambdaA)] <- NA
    phiA[!is.na(phiA)] <- NA
    psiA[!is.na(psiA)] <- NA
    
  } else{
    lambdaA <- lambda
    phiA <- phi
    psiA <- psi
    
  }
  
  lamList <- as.matrix(which(is.na(lambdaA), arr.ind = T))
  il <- nrow(lamList)
  phiList <- as.matrix(which(is.na(phiA), arr.ind = T))
  ip <- il + nrow(phiList)
  psiList <- as.matrix(which(is.na(psiA), arr.ind = T))
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
  vnpsiList[,2] <- paste0(vnpsiList[,1],"~~y", psiList[,2])
  nameList <- rbind(vnlamList, vnphiList, vnpsiList)
  # ========================================================== #
  # ========================================================== #
  # iterate around this function
  fitResults <- matrix(nrow=opt[[2]], ncol=it)
  # progress bar
  if(pb==T) progress_bar <- txtProgressBar(min = 0, max = it, style = 3)
  iter <- 1
  for(iter in 1:it){
    
    # extract iteration information from modList
    x <- modList[iter, ]
    
    # do we need to update lambda?
    if(iter <= il){
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
    if(iter > il & iter <= ip){
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
    if(iter > ip){
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
    
    # use optim() to run simulation
    fit <- optim(inits, get_log_post, control = list(fnscale = -1),
                 hessian = TRUE,
                 sample.data=sampleData, cfa.model=cfaModel)
    param_mean <- fit$par # numerical deriv
    # compute hess at param_mean
    #hess <- numDeriv::hessian(model, param_mean, ...)
    #param_cov_mat <- solve(-hess)
    param_cov_mat <- solve(-fit$hessian)
    
    # scaled covariance matrix (artifically inflate uncertainty)
    scale.cov = opt[[1]]
    A <- diag(scale.cov, nrow=nrow(param_cov_mat), ncol=ncol(param_cov_mat))
    param_cov_mat <- A%*%param_cov_mat%*%t(A)
    
    # sample
    no.samples=opt[[2]]
    fitResults[,iter] <- mcmc(rmvnorm(no.samples, param_mean, param_cov_mat))
    
    if(pb == T) setTxtProgressBar(progress_bar, iter)
  }
  # ========================================================== #
  # ========================================================== #
  
  colnames(fitResults) <- nameList[,2, drop=T]
  
  # Next, standardized (if desired) default
  if(standardize==T){
    # standardize
    obs.var <- extractedLavaan[extractedLavaan[,"mat"]=="theta", ]
    obs.var <- obs.var[which(obs.var$lhs == obs.var$rhs), c("lhs", "est")]
    
    fct.var <- extractedLavaan[extractedLavaan[,"mat"]=="psi", ]
    fct.var <- fct.var[which(fct.var$lhs == fct.var$rhs), c("lhs", "est")]
    
    all.var <- rbind(obs.var, fct.var)
    
    fitResults <- fitResults
    p <- colnames(fitResults)
    i <- 1
    for(i in 1:length(p)){
      unstd <- fitResults[,i]
      
      if(p[i] %like% "=~"){
        pp <- strsplit(p[i], "=~") %>% unlist()
        sigjj <- sqrt(all.var[all.var[,1] == pp[1], 2])
        sigii <- sqrt(all.var[all.var[,1] == pp[2], 2])
        std <- unstd*sqrt(sigjj/sigii) # bollen (1989, p. 349)
      }
      
      if(p[i] %like% "~~"){
        pp <- strsplit(p[i], "~~") %>% unlist()
        sigjj <- sqrt(all.var[all.var[,1] == pp[1], 2])
        sigii <- sqrt(all.var[all.var[,1] == pp[2], 2])
        std <- unstd/(sigjj * sigii) # bollen (1989, p. 349)
      }
      
      fitResults[,i] <- std
    }
  }
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
      pv <- mean(ifelse(abs(x) >= cut.load, 1, 0))
    }
    if(p[i] %like% "~~"){
      pv <- mean(ifelse(abs(x) >= cut.cov, 1, 0))
    }
    sumResults[i, 2] <- pv
  }
  sumResults <- arrange(sumResults, desc(Prob))
  colnames(sumResults) <- c("Parameter","Pr(|theta|>cutoff)", "mean", "sd", "p0.025", "p0.25", "p0.5", "p0.75", "p0.975")
  sumResults[,2:9] <- round(sumResults[,2:9], 3)
  cat("\n")
  if(sum.print==T) print(sumResults, row.names = FALSE)
  
  # convert to data.frame
  fitResults <- as.data.frame(fitResults)
  out <- list(fitResults, sumResults)
  names(out) <- c("All Results", "Summary")
  
  return(out)
}
