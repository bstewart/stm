#a function to evaluate residual dispersion based on
# Taddy's 2012 paper and maptpx code.
checkResiduals <- function(stmobj, documents, tol=.01) {
  beta <- lapply(stmobj$beta$logbeta, exp)
  theta <- stmobj$theta
  index <- stmobj$settings$covariates$betaindex
  K <- stmobj$settings$dim$K
  n <- length(documents)
  phat <- stmobj$settings$dim$V
  d <- n*(K-1) + K*( phat-1 )
  
  doc.resid <- function(doc, theta, beta) {
    q <- theta%*%beta 
    m <- sum(doc[2,])
    Nhat <- sum(q*m > tol)
    x <- rep(0, ncol(beta))
    x[doc[1,]] <- doc[2,]
    out <- sum((x^2 - 2*x*q*m)/(m*q*(1-q))) + sum(m*q/(1-q))
    return(list(out=out, Nhat=Nhat))
  }
  result <- vector(mode="list", length=length(documents))
  Nhat <- 0
  for(i in 1:length(result)) {
    resid <- doc.resid(documents[[i]], theta[i,], beta[[index[i]]])
    result[[i]] <- resid$out
    Nhat <- Nhat + resid$Nhat
  }
  D <- sum(unlist(result))
  df <- Nhat - phat  - d
  sig2 <- D/df
  
  rho <- suppressWarnings(pchisq(D, df=df, lower.tail=FALSE))
  D <- list(dispersion=sig2, pvalue=rho, df=df)
  return(D) 
}
