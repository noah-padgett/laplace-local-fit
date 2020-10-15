# ============================ #
# Bayesian CFa Approximation
#   with Laplace method
# ============================ #
# Created by: R. Noah Padgett
# Created on: 2020-01-10
# 
# Laste Editted: 2020-13
#============================= #
# Adapted from:
# http://www.sumsar.net/blog/2013/11/easy-laplace-approximation/
# by: Rasmus Bååth 
# on: Nov 22nd, 2013
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

X <- lavaan::HolzingerSwineford1939[, paste0('x',1:9)]

nF<-3 # number of factor
p<-ncol(X) # number of variables
N<- nrow(X)# number of individuals

# model specification for factor loading matrix
#   Note: matrix fills column wise
lambdaMod <- matrix(ncol=nF, nrow=p,
                    # x1,x2, x3, x4,x5, x6, x7,x8,x9
                    c(1, 0.55, 0.73, NA,  0,  0, 0,  0, 0,   #f1
                      0,  0,  0, 1, 1.12, 0.93, 0,  0, 0,   #f2
                      0,  0,  0, 0,  0,  0, 1, 1.18, 1.08)) #f3
lambdaMod

# factor covariance matrix (lower diagonal + diagonal)
phiMod <- matrix(nrow=nF,ncol=nF)
diag(phiMod) <- c(0.82, 0.97,  0.38)
phiMod[lower.tri(phiMod)] <-  phiMod[upper.tri(phiMod)] <- c(0.41, 0.262, 0.17);
phiMod

# error variances
psiMod <- diag(nrow=p)
diag(psiMod)<- c(0.55, 1.13, 0.84, 0.37, 0.44, 0.36, 0.80, 0.49, 0.57)
#psiMod[7,1] <- NA
#psiMod[9,1] <- psiMod[7,2] <- NA
psiMod

# intercepts/means
#tauMod <- matrix(ncol=1, nrow=p)

# factor scores
#etaMod <- matrix(ncol=N, nrow=nF)

# store as list
cfaModel <- list(lambdaMod, phiMod, psiMod) #, tauMod, etaMod

# starting values
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
sv[1:(lam.num)] <- rep(0.25,lam.num)
  if(lam.num==0){ sv.n1 <- NA }else{ sv.n1 <-   paste0('lambda', 1:lam.num)}
sv[(lam.num+1):(lam.num+dphi.num)]<- runif(dphi.num, 0.05, 1)
  if(dphi.num==0){ sv.n2 <- NA }else{ sv.n2 <-   paste0('dphi', 1:dphi.num)}
sv[(lam.num+dphi.num+1):(lam.num+phi.num)]<- runif(odphi.num, -.1, 0.1)
  if(odphi.num==0){sv.n3 <- NA }else{ sv.n3 <-  paste0('odphi', 1:odphi.num)}
sv[(lam.num+phi.num+1):(lam.num+phi.num+dpsi.num)] <- rep(0.2, dpsi.num)
  if(dpsi.num==0){ sv.n4 <- NA }else{ sv.n4 <-  paste0('dpsi', 1:dpsi.num)}
sv[(lam.num+(phi.num + dpsi.num)+1):(lam.num+phi.num+psi.num)] <- runif(odpsi.num, -.05, 0.05)
  if(odpsi.num==0){ sv.n5 <- NA }else{ sv.n5 <-  paste0('odpsi', 1:odpsi.num)}
#sv[(lam.num+phi.num+psi.num+1):(lam.num+phi.num+psi.num+tau.num)] <-sample(unlist(X), tau.num, replace=T)
#sv[(lam.num+phi.num+psi.num+tau.num+1):(lam.num+phi.num+psi.num+tau.num+eta.num)] <-rnorm(eta.num, 0, 1)
  names(sv) <- na.omit(c(sv.n1, sv.n2, sv.n3, sv.n4, sv.n5)) #, paste0('tau', 1:tau.num), paste0('eta', 1:eta.num))

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

f <- function(m) {
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
  phi <- f(phi) # map lower to upper
  
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
  psi <- f(psi)
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


repeat_col<-function(x,n){
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


inits <- sv # initial values
names(inits)<- names(sv)

samples <- laplace_approx(model, inits, 10000, X=X, cfaModel=cfaModel)
coda::densplot(samples)
summary(samples)


HS.model <- ' 
visual  =~ x1 + x2 + x3 + x4
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
'

fit <- cfa(HS.model, data = HolzingerSwineford1939)
summary(fit, fit.measures = TRUE)
lavResiduals(fit, type='cor')



