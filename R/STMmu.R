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
    gamma <- t(do.call(rbind,gamma))
    mu<- t(covar%*%gamma)
    return(list(mu=mu, gamma=gamma))
  }
  
  #Gamma Lasso
  if(mode=="L1") {
    out <- glmnet(x=covar[,-1], y=lambda, family="mgaussian", alpha=enet)
    unpack <- unpack.glmnet(out, nobs=nrow(covar), ic.k=2)
    gamma <- rbind(unpack$intercept, unpack$coef)
    mu <- t(covar%*%gamma)
    return(list(mu=mu, gamma=gamma))
  }
}

#Variational Linear Regression with a Gamma hyperprior Defaults to a broad (1,1)
# (See Drugowitsch 2013 or Murphy)
vb.variational.reg <- function(Y,X, a0=1, b0=1) {
  X <- X[,-1,drop=FALSE]
  Xm <- colMeans(X)
  Ym <- mean(Y)
  X.uncentered <- X
  Y.uncentered <- Y
  X <- sweep(X,2, STATS=Xm,FUN="-")
  Y <- Y - Ym
  
  c0 <- 1
  d0 <- 1
  Xcorr <- crossprod(X)
  XYcorr <- crossprod(X,Y) 
  an <- a0 + nrow(X)/2
  cn <- c0 + ncol(X)/2
  D <- ncol(X); N <- nrow(X)
  Ea <- c0/d0
  w <- rep(0, ncol(X))
  converge <- 1000
  while(converge>.0001) {
    w.old <- w
    invV <- as.numeric(Ea)*diag(1, nrow=D) + Xcorr
    V <- solve(invV)
    w <- V%*%XYcorr
    # parameters of noise model (an remains constant)
    sse <- sum((X %*% w - Y)^ 2)
    bn <- b0 + 0.5 * (sse + Ea %*%crossprod(w))
    Et <- an / bn
    #hyperparameters of covariance prior (cn remains constant)
    dn <- d0 + 0.5 * (Et %*% crossprod(w) + sum(diag(V)))
    Ea = cn / dn
    converge <- sum(abs(w-w.old))
  }
  w0 <- mean(Y.uncentered) - mean(X.uncentered%*%w)
  out <- c(w0,w)
  return(out)
}
