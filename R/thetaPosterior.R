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

#####
#Recompute the local covariance matrix
# then draw from the posterior distribution
##
# NB: this is slightly complicated by the need to take Newton Steps
thetapost.local <- function(model, documents, nsims) {
  #define some parameters
  sigma <- model$sigma
  siginv <- model$invsigma
  mu <- model$mu$mu
  lambda <- model$eta
  beta <- lapply(model$beta$logbeta, exp)
  betaindex <- model$settings$covariates$betaindex
    
  out <- vector(mode="list",length=nrow(lambda)) 
  for(i in 1:nrow(lambda)) {
    #for each document we take the following steps:
    # 1) check if inverse is a legitimate cov matrix
    # 2) if not, tighten with optim and check again
    # 3) if still not, do newton steps
    doc.ct <- documents[[i]][2,]
    eta <- lambda[i,]
    theta <- model$theta[i,]
    #get just the necessary columns of beta
    doc.beta <-   beta[[betaindex[i]]][,documents[[i]][1,]]
    hess <- ln.hess(eta, theta, doc.beta, doc.ct, siginv)    
    nu <- try(chol2inv(chol.default(hess)),silent=TRUE)
    if(class(nu)=="try-error") {
      # if it failed we try tightening by taking a BFGS step
      optim.out <- optim(par=eta, fn=lhoodcpp, gr=gradcpp,
                         method="BFGS", control=list(maxit=500),
                         doc.ct=doc.ct, mu=mu[,i],
                         siginv=siginv, beta=doc.beta, Ndoc=sum(doc.ct))
      eta <- optim.out$par
      theta <- softmax(c(eta,0))
      hess <- ln.hess(eta, theta, doc.beta, doc.ct, siginv) 
      nu <- try(chol2inv(chol.default(hess)),silent=TRUE)
      if(class(nu)=="try-error") {
        #oh man, it failed again?  Well now we try some really expensive newton optimization
        newton.out <- newton(eta, doc.ct, mu[,i], siginv, doc.beta, hess, max.its=1000)
        #if the newton even failed it forces an answer using nearPD so we do have a solution here.
        nu <- newton.out$nu
        eta <- newton.out$eta
      }
    }
    #take the draws
    mat <- rmvnorm(nsims, eta,nu)
    #add the 0
    mat <- cbind(mat,0)
    #project to the simplex and save
    out[[i]] <- exp(mat - row.lse(mat))
  }
  return(out)
}

#Full Newton Optimization
# a costly Newton optimization search strategy.
# this is used when standard Quasi-Newton methods
# have failed.
newton <- function(eta, doc.ct, mu, siginv, beta, 
                   hess, max.its=1000) {
  its <- 0
  search <- function(x, dir, eta, ...) {
    lhoodcpp(eta=(eta + x*dir), ...)
  }
  while(its < max.its) {
    #compute the search direction
    dir <- -solve(hess)%*%gradcpp(eta, doc.ct, mu, siginv, beta)
    #line search
    maxint <- 2
    opt <- optimize(f=search, interval=c(-2,maxint), dir=dir, 
                    eta=eta,doc.ct=doc.ct, mu=mu,
                    siginv=siginv, beta=beta, maximum=FALSE)
    while(opt$objective > lhoodcpp(eta, doc.ct, mu, siginv, beta)) {
      #re-assess at a smaller point
      maxint <- min(maxint/2, opt$minimum - .00001*opt$minimum)
      opt <- optimize(f=search, interval=c(0,maxint), dir=dir, 
                      eta=eta,doc.ct=doc.ct, mu=mu,
                      siginv=siginv, beta=beta, maximum=FALSE)
    }
    
    eta <- eta + opt$minimum*dir
    expeta <- c(exp(eta),1)
    theta <- expeta/sum(expeta)
    hess <- ln.hess(eta, theta, beta, doc.ct, siginv) 
    nu <- try(chol2inv(chol.default(hess)), silent=TRUE)
    if(class(nu)!="try-error") {
      break
    } else {
      its <- its + 1
    }
  }
  if(its==max.its) {
    warning(sprintf("Document covariance matrix didn't converge after %i Newton iterations. Using nearest positive definite matrix to inverse hessian.", max.its))
    nu <- nearPD(solve(hess))
  }
  return(list(eta=eta,nu=nu))
}



ln.hess <- function(eta, theta, beta,doc.ct, siginv) {
  expeta <- c(exp(eta),1)
  EB <- expeta*beta #calculate exp(eta)\beta for each word
  EB <- t(EB)/colSums(EB) #transpose and norm by (now) the row
  phi <- EB*(doc.ct) #multiply through by word count
  phisums <- colSums(phi)
  EB <- EB*sqrt(doc.ct) #set up matrix to take the cross product
  Ndoc <- sum(doc.ct)
  hess <- -((diag(phisums) - crossprod(EB)) - 
              Ndoc*(diag(theta) - theta%o%theta))[1:length(eta),1:length(eta)] + siginv
  return(hess)
}
