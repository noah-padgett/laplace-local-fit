# ============================ #
# Bayesian CFA Approximation
#   with Laplace method
#   Generalized version
# ============================ #
# Created by: R. Noah Padgett
# Created on: 2020-01-10
# 
# Laste Editted: 2020-02-18
# ============================ #


library(coda)
library(mvtnorm)
library(lavaan)
library(data.table)

# Laplace approximation
laplace_approx <- function(model, inits, no_samples, ...) {
  fit <- optim(inits, model, control = list(fnscale = -1), hessian = TRUE,
               ...)
  param_mean <- fit$par
  param_cov_mat <- solve(-fit$hessian)
  mcmc(rmvnorm(no_samples, param_mean, param_cov_mat))
}


HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939)
summary(fit, fit.measures = TRUE)

# extract model from lavaan
matrixConversion <- function(x, y, est) {
  Z<- matrix(NA, nrow=max(x), ncol=max(y))
  for(i in 1:length(est)){
    Z[x[i], y[i]] <- est[i]
  }
  Z
}



extractedLavaan <- lavMatrixRepresentation(partable(fit))

factNames <- unique(extractedLavaan[extractedLavaan[,"mat"]=="lambda", "lhs"])
varNames <- unique(extractedLavaan[extractedLavaan[,"mat"]=="lambda", "rhs"])
# extract factor loading matrix
lambda <- extractedLavaan[ extractedLavaan$mat == "lambda" ,]
lambda <- matrixConversion(lambda$row, lambda$col, lambda$est)
colnames(lambda) <- factNames
rownames(lambda) <- varNames
# extract factor covariance matrix
phi <- extractedLavaan[ extractedLavaan$mat == "psi" ,]
phi <- matrixConversion(phi[,'row'], phi[,'col'], phi[,'est'])
phi <- up2full(phi)
colnames(phi) <- rownames(phi) <- factNames
# extract error covariance matrix
psi <- extractedLavaan[ extractedLavaan$mat == "theta" ,]
psi <- matrixConversion(psi[,'row'], psi[,'col'], psi[,'est'])
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




# set up first variable to apply method to
i <- 1
x <- modList[i,]
lambdaMod <- lambda
phiMod <- phi
psiMod <- psi
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

# starting values
get_starting_values <- function(model){
  lambdaMod <- model[[1]]
  phiMod <- model[[2]]
  psiMod <- model[[3]]

  # get length of each model element
  lam.num <- sum(is.na(c(lambdaMod))==T)
  phi.num <- sum(is.na(c(phiMod))==T)
  dphi.num <- sum(is.na(diag(phiMod))==T)
  odphi.num <- sum(is.na(phiMod[lower.tri(phiMod)])==T)
  psi.num <- sum(is.na(c(psiMod))==T)
  dpsi.num <- sum(is.na(diag(psiMod))==T)
  odpsi.num <- sum(is.na(psiMod[lower.tri(psiMod)])==T)
  # tau.num <- p
  # eta.num <- length(etaMod)
  
  k<-lam.num+phi.num+psi.num#+tau.num+eta.num
  sv<-numeric(k)
  # generate starting values
  sv.n1 <- sv.n2 <- sv.n3 <- sv.n4 <- sv.n5 <- NA
  if(lam.num==0){ 
    sv.n1 <- NA 
  }else{
    sv[1:(lam.num)] <- runif(lam.num, 0.6, 0.8)
    sv.n1 <-   paste0('lambda', 1:lam.num)
  }
  if(dphi.num==0){
    sv.n2 <- NA 
  }else{
    sv[(lam.num+1):(lam.num+dphi.num)]<- runif(dphi.num, 0.05, 1)
    sv.n2 <-   paste0('dphi', 1:dphi.num)
  }
  if(odphi.num==0){
    sv.n3 <- NA 
  }else{ 
    sv[(lam.num+dphi.num+1):(lam.num+phi.num)]<- runif(odphi.num, -.1, 0.1)
    sv.n3 <-  paste0('odphi', 1:odphi.num)
  }
  
  if(dpsi.num==0){
    sv.n4 <- NA 
  }else{ 
    sv[(lam.num+phi.num+1):(lam.num+phi.num+dpsi.num)] <- rep(0.2, dpsi.num)
    sv.n4 <-  paste0('dpsi', 1:dpsi.num)
  }
  
  if(odpsi.num==0){
    sv.n5 <- NA
  }else{ 
    sv[(lam.num+(phi.num + dpsi.num)+1):(lam.num+phi.num+psi.num)] <- runif(odpsi.num, -.05, 0.05)
    sv.n5 <-  paste0('odpsi', 1:odpsi.num)
  }
  #sv[(lam.num+phi.num+psi.num+1):(lam.num+phi.num+psi.num+tau.num)] <-sample(unlist(X), tau.num, replace=T)
  #sv[(lam.num+phi.num+psi.num+tau.num+1):(lam.num+phi.num+psi.num+tau.num+eta.num)] <-rnorm(eta.num, 0, 1)
  names(sv) <- na.omit(c(sv.n1, sv.n2, sv.n3, sv.n4, sv.n5)) #, paste0('tau', 1:tau.num), paste0('eta', 1:eta.num))
  return(sv)
}
# trace function
trace <- function(A) {
  n <- dim(A)[1] # get dimension of matrix
  tr <- 0 # initialize trace value
  
  # Loop over the diagonal elements of the supplied matrix and add the element to tr
  for (k in 1:n) {
    l <- A[k,k]
    tr <- tr + l
  }
  return(tr[[1]])
}
# or one could do sum(diag(A))

up2full <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}
low2full <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

# need to use function to easily manipulate CFA model specifications

use_cfa_model <- function(theta, X, model){
  # Compue sample statistics
  p<-ncol(X)
  S<-cov(X)
  
  # unpack model
  lambda <- model[[1]]
  phi <- model[[2]]
  psi <- model[[3]]
  #tau <- model[[4]]
  #eta <- model[[5]]
  
  # number factor loadings
  lam.num <- length(which(is.na(lambda)))
  lambda[which(is.na(lambda))] <- theta[1:lam.num]
  nF = ncol(lambda)
  # number elements in factor (co)variance matrix
  phi.num <- length(which(is.na(phi)))
  dphi.num <- sum(is.na(diag(phiMod))==T)
  odphi.num <- sum(is.na(phiMod[lower.tri(phiMod)])==T)
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



repeat_col <- function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

get_prior_dens <- function(pvalue, pname){
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


model <- function(p, X, cfaModel) {
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



# running the method
inits <- get_starting_values(cfaModel) # initial values
X <- lavaan::HolzingerSwineford1939[, paste0('x',1:9)]


samples <- laplace_approx(model, inits, 1000, X=X, cfaModel=cfaModel)
coda::densplot(samples)
summary(samples)


