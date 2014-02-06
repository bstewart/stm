#Optimization for the Global Covariance Matrix
opt.sigma <- function(nu, lambda, mu, sigprior) {  

  #find the covariance
  if(ncol(mu)==1) {
    covariance <- crossprod(sweep(lambda, 2, STATS=as.numeric(mu), FUN="-"))
  } else {
    covariance <- crossprod(lambda-t(mu)) 
  }
  sigma <- (covariance + nu)/nrow(lambda) #add to estimation variance
  sigma <- diag(diag(sigma),nrow=nrow(nu))*sigprior + (1-sigprior)*sigma #weight by the prior
  return(sigma)
}


