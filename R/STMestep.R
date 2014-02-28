#Local Logistic-Normal Inference Using Standard Batch Scheduling

#Input: Documents and Key Global Parameters
#Output: Sufficient Statistics

# Approach:
# First we pre-allocate memory, and precalculate where possible.
# Then for each document we:
#  (1) get document-specific priors, 
#  (2) infer doc parameters, 
#  (3) update global sufficient statistics
# Then the sufficient statistics are returned.

estep.LN <- function(documents, beta.index, logbeta, mu, sigma, verbose, lambdacurrent=NULL) {
  #quickly define useful constants
  V <- ncol(logbeta[[1]])
  K <- nrow(logbeta[[1]])
  A <- length(unique(beta.index))
  N <- length(documents)
  if(verbose) {
    if(N>100) {
      ctevery <- floor(N/100)
    } else {
      ctevery <- 1
    }
  }
  # 1) Initialize Sufficient Statistics 
  sigma.ss <- diag(0, nrow=(K-1))
  beta.ss <- vector(mode="list", length=A)
  for(i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
    
  priors <- list(mu=NULL,
                 sigmaentropy=(.5*determinant(sigma, logarithm=TRUE)$modulus[1]),
                 siginv=solve(sigma), 
                 logbeta=NULL)
  if(ncol(mu)==1) {
    updateMu <- FALSE
    priors$mu <- as.vector(mu) #converting to a vector here and below solves dimensionality issues in matrix mult.
  } else {
    updateMu <- TRUE
  }
  
  # 2) Document Scheduling
  for(i in 1:N) {
    # a) bundle parameters
    doc <- documents[[i]]
    words <- doc[1,]
    aspect <- beta.index[i]
    if(updateMu) priors$mu <- as.vector(mu[,i])
    priors$logbeta <- logbeta[[aspect]][,words,drop=FALSE]
    
    if(is.null(lambdacurrent)) {
      init <- priors$mu
    } else {
      init <- lambdacurrent[i,]
    }
    # b) infer local latent variables 
    doc.results <- inferdoc.logisticnormal(doc,priors,init)
    
    # c) update sufficient statistics 
    sigma.ss <- sigma.ss + doc.results$eta$nu
    beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
    bound[i] <- doc.results$bound
    lambda[[i]] <- doc.results$eta$lambda
    if(verbose && i%%ctevery==0) cat(".")
  }
  if(verbose) cat("\n") #add a line break for the next message.
  
  #3) Combine and Return Sufficient Statistics
  lambda <- do.call(rbind, lambda)
  return(list(sigma=sigma.ss, beta=beta.ss, bound=bound, lambda=lambda))
}
