## Optimization for Global Parameters over Doc-Topic Proportions

#main method up top, regression-implementations below.
opt.mu <- function(lambda, mode=c("CTM","Pooled", "L1"), covar=NULL, enet=NULL) {

  #When there are no covariates we use the CTM method
  if(mode=="CTM") {
    mu <- matrix(colMeans(lambda), ncol=1)
    return(list(mu=mu, gamma=NULL))
  }
  
  #Variational Linear Regression with a Gamma hyperprior
  if(mode=="Pooled") {
    gamma <- vector(mode="list",length=ncol(lambda))
    for (i in 1:ncol(lambda)) {
      gamma[[i]] <- vb.variational.reg(Y=lambda[,i], X=covar) 
    }
    gamma <- do.call(cbind,gamma)
    mu<- t(covar%*%gamma)
    #if its not a regular matrix,coerce it as it won't be sparse.
    if(!is.matrix(mu)) {
      mu <- as.matrix(mu)
    }
    return(list(mu=mu, gamma=gamma))
  }
  
  #Lasso
  if(mode=="L1") {
    out <- glmnet(x=covar[,-1], y=lambda, family="mgaussian", alpha=enet)
    unpack <- unpack.glmnet(out, nobs=nrow(covar), ic.k=2)
    gamma <- rbind(unpack$intercept, unpack$coef)
    mu <- t(covar%*%gamma)
    if(!is.matrix(mu)) {
      mu <- as.matrix(mu)
    }
    return(list(mu=mu, gamma=gamma))
  }
}

#Variational Linear Regression with a Half-Cauchy hyperprior 
# (Implementation based off the various LMM examples from Matt Wand)
# This code is intended to be passed a Matrix object
vb.variational.reg <- function(Y,X, b0=1, d0=1) {
  Xcorr <- crossprod(X)
  XYcorr <- crossprod(X,Y) 
  
  an <- (1 + nrow(X))/2
  D <- ncol(X)
  N <- nrow(X)
  w <- rep(0, ncol(X))
  error.prec <- 1 #expectation of the error precision
  converge <- 1000
  cn <- ncol(X) # - 1 for the intercept and +1 in the update cancel
  dn <- 1
  Ea <- cn/dn #expectation of the precision on the weights
  ba <- 1
  
  while(converge>.0001) {
    w.old <- w
    
    #add the coefficient prior.  Form depends on whether X is a Matrix object or a regular matrix.
    if(is.matrix(X)) {
      ppmat <- diag(x=c(0, rep(as.numeric(Ea), (D-1))),nrow=D) 
    } else {
      ppmat <- Diagonal(n=D, x=c(0, rep(as.numeric(Ea), (D-1))))
    }
    invV <- error.prec*Xcorr + ppmat
    #if its a plain matrix its faster to use the cholesky, otherwise just use solve
    if(is.matrix(invV)) {
      V <- chol2inv(chol(invV))
    } else {
      #Matrix package makes this faster even when its non-sparse
      V <- solve(invV)     
    }
    w <- error.prec*V%*%XYcorr
    
    # parameters of noise model (an remains constant)
    sse <- sum((X %*% w - Y)^ 2)
    bn <- .5*(sse + sum(diag(Xcorr%*%V))) + ba
    error.prec <- an/bn
    ba <- 1/(error.prec + b0)
    
    #subtract off the intercept while working out the hyperparameters
    # for the coefficients
    w0 <- w[1]
    w <- w[-1]
    da <- 2/(Ea + d0)
    dn <- 2*da + (crossprod(w) + sum(diag(V)[-1]))
    Ea <- cn / dn
    #now combine the intercept back in 
    w <- c(w0,w)
    converge <- sum(abs(w-w.old))
  }
  return(w)
}


