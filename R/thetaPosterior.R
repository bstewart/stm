#The exported function with workhorse functions below
thetaPosterior <- function(model, nsims=100, type=c("Global", "Local"), documents=NULL) {
  type <- match.arg(type)
  if(type=="Local" & is.null(documents)) stop("Documents must be provided to perform local theta uncertainty calculations.")
  switch(type,
         Global=thetapost.global(model, nsims),
         Local=thetapost.local(model, documents, nsims))
}
  
#Take nsims draws from the variational posterior over theta using 
# a global approximation to the covariance matrix
thetapost.global <- function(model, nsims) {
  Sigma <- model$sigma
  mu <- model$mu$mu
  lambda <- model$eta
  
  ##
  #Calculate approximate global covariance matrix
  ##
  #the general idea is to subtract off the component of the Sigma update arising from deviation to the
  # global mean leaving us with the component that comes from nu
  if(ncol(mu)==1) { 
    #if there is only one global mean vector we avoid storing them all, thus calc done with a sweep
    covariance <- crossprod(sweep(lambda, 2, STATS=as.numeric(mu), FUN="-"))
  } else {
    #the typical calculation when there are frequency covariates
    covariance <- crossprod(lambda-t(mu)) 
  }
  #rescale by the number of documents
  covariance <- covariance/nrow(lambda)
  #subtract off the effect from the global covariance
  Sigma <- Sigma - covariance
  
  ##
  #Sample and Project to Simplex
  ##
  choleskydecomp <- chol(Sigma)
  out <- vector(mode="list",length=nrow(lambda)) 
  for (i in 1:length(out)) {
    mat <- rmvnorm(nsims, lambda[i,],Sigma,choleskydecomp)
    mat <- cbind(mat, 0)
    out[[i]] <- exp(mat - row.lse(mat))
  }
  return(out)
}
  
thetapost.local <- function(model, documents, nsims) {
  #the local approximation
  Sigma <- model$sigma
  siginv <- model$invsigma
  mu <- model$mu$mu
  lambda <- model$eta
  logbeta <- model$beta$logbeta
  betaindex <- model$settings$covariates$betaindex
    
    #Note to ensure positive definiteness we essentially have to do another E-step
    calc.nu <- function(init, doc.ct, doc.mu, doc.logbeta) {
      optim.out <- optim(par=init, fn=eta.likelihoodjoint, gr=eta.gradientjoint,
                         method="L-BFGS-B",hessian=FALSE,
                         doc.ct=doc.ct, mu=doc.mu,
                         siginv=siginv, logbeta=doc.logbeta)
	  hessian <- eta.hessianjoint(optim.out$par, 
					     doc.ct=doc.ct, mu=doc.mu,
					     siginv=siginv, logbeta=doc.logbeta, 
						 Ndoc=sum(doc.ct))
      est <- optim.out$par
      nu <- solve(hessian$hessian)
      return(list(est=est, nu=nu))
    }
    
    out <- vector(mode="list",length=nrow(lambda)) 
    for(i in 1:nrow(lambda)) {
      init <- lambda[i,]
      
      if(ncol(mu)==1) {
        doc.mu <- as.numeric(mu)
      } else {
        doc.mu <- as.numeric(mu[,i])
      }
      doc <- documents[[i]]
      doc.logbeta <- logbeta[[betaindex[i]]][,doc[1,],drop=FALSE]
      doc.ct <- doc[2,drop=FALSE]
      params <- calc.nu(init, doc.ct, doc.mu, doc.logbeta)
      mat <- rmvnorm(nsims, params$est,params$nu)
      mat <- cbind(mat,0)
      out[[i]] <- exp(mat - row.lse(mat))
    }
    return(out)
  }