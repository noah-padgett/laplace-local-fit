# ============================ #
# Bayesian CFA Approximation
#   with Laplace method
#   
#   Utility Functions
#
# ============================ #
# Created by: R. Noah Padgett
# Created on: 2020-02-18
# 
# Laste Editted: 2020-02-18
# ============================ #

# ========================================== #
# ========================================== #
#   function: convert2matrix()
# ========================================== #
# use: takes vectors of positions in matrix
#       and values to create a possible 
#       sparse matrix
#
# arguments:
# theta - vector of parameters being optimized
# X     - obversed data
# model - list of model components (lambda...)
#
convert2matrix <- function(x, y, est) {
  Z<- matrix(NA, nrow=max(x), ncol=max(y))
  for(i in 1:length(est)){
    Z[x[i], y[i]] <- est[i]
  }
  Z
}

# ========================================== #
# ========================================== #
#   function: get_starting_values()
# ========================================== #
# use: obtain starting values, will work
#       for any parameters
#
# arguments:
# model - list of model components (lambda...)
#
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
# ========================================== #
# ========================================== #
#   function: trace()
# ========================================== #
# use: compute trace of matrix
#
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
# ========================================== #
# ========================================== #
#   function: XX2full()
# ========================================== #
# use: convert lower or upper matrix to full
up2full <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}
low2full <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

# ========================================== #
# ========================================== #
#   function: repeat_cols()
# ========================================== #
# use: take in vector and converts into long matrix
repeat_col <- function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


# ========================================== #
# ========================================== #
#   function: progress_bar()
# ========================================== #
# use: makes the progress bar 
progress_bar <- txtProgressBar(min = 0, max = 2, style = 3)
