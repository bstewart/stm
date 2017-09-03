# The exported function with workhorse functions below

#' Draw from Theta Posterior
#' 
#' Take random draws from the variational posterior for the document-topic
#' proportions. This is underlying methodology for \code{\link{estimateEffect}}
#' 
#' This function allows the user to draw samples from the variational posterior
#' distribution over the normalized document-topic proportions, theta. The
#' function \code{\link{estimateEffect}} provides a user-friendly interface for
#' running regressions using samples from the posterior distribution.  When the
#' user wants to do something not covered by that function, the code here
#' provides easy access to uncertainty in the model.
#' 
#' In order to simulate from the variational posterior for theta we take draws
#' from the variational distribution for eta (the unnormalized topic
#' proportions) and then map them to the simplex.  Each document in the corpus
#' has its own mean vector (eta) and covariance matrix (nu).  Because the
#' covariance matrices can be large we do not store them in the model objects.
#' We offer two approximations to the covariance matrix: Global and Local.  The
#' Global method constructs a single approximate covariance matrix which is
#' then shared by all documents.  This approach is very fast and does not
#' require access to the original documents.  For highly aggregated quantities
#' of interest this often produces similar results to the Local method.
#' 
#' The Local method steps through each document in sequence and calculates the
#' covariance matrix.  If the model has not converged, this matrix can be
#' undefined and so we perform document level inference until the estimate
#' stabilizes.  This means that under the Local method both the covariance and
#' the mean of the variational distribution are recalculated.  It also means
#' that calling this option with Local specified will take approximately as
#' long as a standard E-step of \code{\link{stm}} for the same data and
#' possibly longer.  Because the memory requirements would be extreme for large
#' K, we calculate one document at a time, discarding the covariance matrix
#' before proceeding to the next document.  Thus, if your computer has
#' sufficient memory it is dramatically more computationally efficient to draw
#' all the samples you may require at once rather than taking one sample at a
#' time.
#' 
#' The output for both methods is a list with number of elements equal to the
#' number of documents. Each element is a matrix with nsims rows and K columns.
#' Be careful to ensure that you have sufficient memory before attempting this
#' with a large number of simulations, documents or topics.
#' 
#' @param model An \code{STM} object created by \code{\link{stm}}
#' @param nsims The number of draws from the variational posterior.  See
#' details below.
#' @param type A choice of two methods for constructing the covariance
#' approximation the \code{"Global"} approximation and the \code{"Local"}
#' approximation.  See details below.
#' @param documents If \code{type="Local"}, the documents object used in the
#' original \code{\link{stm}} call should be passed here.
#' @seealso \code{\link{estimateEffect}}
#' @examples
#' #global approximation
#' draws <- thetaPosterior(gadarianFit, nsims = 100)
#' @export
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
    doc.beta <-   beta[[betaindex[i]]][,documents[[i]][1,], drop=FALSE]
    hess <- ln.hess(eta, theta, doc.beta, doc.ct, siginv)    
    nu <- try(chol2inv(chol.default(hess)),silent=TRUE)
    if(class(nu)=="try-error") {
    
      # if it failed we try tightening by taking a BFGS step
      if (NCOL(mu) == 1) mu.i <- mu else mu.i <- mu[,i]
      optim.out <- optim(par = eta, fn = lhoodcpp, gr = gradcpp, 
                         method = "BFGS", control = list(maxit = 500), 
                         doc_ct = doc.ct, mu = mu.i, siginv = siginv, 
                         beta = doc.beta)

      eta <- optim.out$par
      theta <- softmax(c(eta,0))
      hess <- ln.hess(eta, theta, doc.beta, doc.ct, siginv) 
      nu <- try(chol2inv(chol.default(hess)),silent=TRUE)
      if(class(nu)=="try-error") {
      
        #oh man, it failed again?  Well now we try some really expensive newton optimization
        newton.out <- newton(eta, doc.ct, mu.i, siginv, doc.beta, hess, max.its=1000)
        
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
    dir <- -solve(hess)%*%gradcpp(eta=eta, doc_ct=doc.ct, 
                                  mu=mu, siginv=siginv, beta=beta)

    #line search
    maxint <- 2
    opt <- optimize(f=search, interval=c(-2,maxint), dir=dir, 
                    eta=eta,doc_ct=doc.ct, mu=mu,
                    siginv=siginv, beta=beta, maximum=FALSE)
    while(opt$objective > lhoodcpp(eta, beta, doc.ct, mu, siginv)) {
      #re-assess at a smaller point
      maxint <- min(maxint/2, opt$minimum - .00001*opt$minimum)
      opt <- optimize(f=search, interval=c(0,maxint), dir=dir, 
                      eta=eta,doc_ct=doc.ct, mu=mu,
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