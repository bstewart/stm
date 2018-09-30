# Serial Implementation of the E-Step
estepSerial <- function(N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv, verbose) {
  
  sigma.ss <- diag(0, nrow=(K-1))
  beta.ss <- vector(mode="list", length=A)
  for(i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
  ctevery <- ifelse(N>100, floor(N/100), 1)
  
  for (i in 1:N) {
    doc <- documents[[i]]
    words <- doc[1,]
    aspect <- beta.index[i]
    init <- lambda.old[i,]
    if (update.mu) {
      mu.i <- mu[, i]
    } else {
      mu.i <- as.numeric(mu)
    }
    beta.i <- beta[[aspect]][,words,drop=FALSE]
    
    #infer the document
    doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, doc=doc, sigmaentropy=sigmaentropy)
    
    # update sufficient statistics
    sigma.ss <- sigma.ss + doc.results$eta$nu
    beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
    bound[i] <- doc.results$bound
    lambda[[i]] <- c(doc.results$eta$lambda)
    if(verbose && i%%ctevery==0) cat(".")
  }
  if(verbose) cat("\n") #add a line break for the next message.
  
  lambda <- do.call(rbind, lambda)
  list(sigma=sigma.ss, beta=beta.ss, bound=bound, lambda=lambda)
}


# Parallel Implementation of the E-Step
estepParallel <- function(N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv, verbose, cores=1) {
  
  # INITIALIZATION OF COMBINED RESULTS
  beta.ss <- vector(mode="list", length=A)
  for (i in 1:A) beta.ss[[i]] <- matrix(0, nrow=K, ncol=V)
  init <- list(
    sigma.ss = diag(0, nrow=(K-1)),
    beta.ss = beta.ss,
    bound = vector(length=N),
    lambda = vector("list", length=N)
  )

  # The combine function in which we accumulate results as they're available
  # R represents the combined result object, r the individual result object
  combineFn <- function(R, r) {
    # update sufficient statistics
    R$sigma.ss <- R$sigma.ss + r$sigma.ss
    for (i in length(R$beta.ss)) {
      R$beta.ss[[i]] =  R$beta.ss[[i]] + r$beta.ss[[i]]
    }
    R$bound[r$doc.ids] <- r$bound[r$doc.ids]
    R$lambda[r$doc.ids] <- r$lambda[r$doc.ids]
	
	  # Note: The combine function needs to return the accumulated result object
    R
  }
  
  if (verbose) cat("Starting Parallel E-Step\n")
  
  # Create groups of doc.ids that we will schedule to the parallel workers
  doc.id.groups <- base::split(seq_len(N), sample(rep(seq_len(cores), length=N)))
  
  # We cannot seem to use the foreach::`%dopar%` syntax directly without capturing the operator locally first
  `%dopar%` <- foreach::`%dopar%`
  res <- foreach::foreach (doc.ids = doc.id.groups, .combine = combineFn, .multicombine = FALSE, .init = init) %dopar% {
    estepParallelBlock(doc.ids, N, K, A, V, documents[doc.ids], beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv)
  }
  
  lambda <- do.call(rbind, res$lambda)
  list(sigma=res$sigma.ss, beta=res$beta.ss, bound=res$bound, lambda=lambda)
}

# The E-Step code block that can be scheduled in parallel to other similar blocks
# Each estepParallelBlock is passed a non-overlapping subset of the 'documents' matrix,
# corresponding to doc.ids, which are used to index into other data structures that are also passed in
# The doc.ids are also part of the return values of this block so that the combine function knows what to do.
estepParallelBlock <- function(doc.ids, N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv) {
  
  sigma.ss <- diag(0, nrow=K-1)
  beta.ss <- vector(mode='list', length=A)
  for(i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow=K, ncol=V)
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)

  if(!update.mu) mu.i <- as.numeric(mu)
  
  for (i in 1:length(doc.ids)) {
    doc <- documents[[i]]
	  doc.id <- doc.ids[i]
    words <- doc[1,]
    aspect <- beta.index[doc.id]
    init <- lambda.old[doc.id,]
    if (update.mu) mu.i <- mu[, doc.id]
    beta.i <- beta[[aspect]][, words, drop=FALSE]
    
    #infer the document
    doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, doc=doc, sigmaentropy=sigmaentropy)
    
    # update sufficient statistics
    sigma.ss <- sigma.ss + doc.results$eta$nu
    beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
    bound[doc.id] <- doc.results$bound
    lambda[[doc.id]] <- c(doc.results$eta$lambda)
  }
  list(doc.ids=doc.ids, sigma.ss=sigma.ss, beta.ss=beta.ss, bound=bound, lambda=lambda)
}

#E-Step for a Document Block
#[a relatively straightforward rewrite of previous
# code with a focus on avoiding unnecessary computation.]

#Input: Documents and Key Global Parameters
#Output: Sufficient Statistics

# Approach:
# First we pre-allocate memory, and precalculate where possible.
# Then for each document we:
#  (1) get document-specific priors, 
#  (2) infer doc parameters, 
#  (3) update global sufficient statistics
# Then the sufficient statistics are returned.

#Let's start by assuming its one beta and we may have arbitrarily subset the number of docs.
estep <- function(documents, beta.index, update.mu, #null allows for intercept only model  
                       beta, lambda.old, mu, sigma, 
                       verbose, cores=1, sigma.round=NULL, beta.round=NULL) {
  
  #quickly define useful constants
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  if (!update.mu) mu.i <- as.numeric(mu)
  
  # Precalculate common components
  sigobj <- try(chol.default(sigma), silent=TRUE)
  if (class(sigobj)=="try-error") {
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
    sigmaentropy <- sum(log(diag(sigobj)))
    siginv <- chol2inv(sigobj)
  }
  
  if (cores>1) {
    results <- estepParallel(N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv, verbose, cores)
  } else {
    results <- estepSerial(N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv, verbose)
  }
  
  # E-estep may yield slightly different results depending on whether it's run in serial or parallel,
  # and how the work is split between workers. This is because of the inherent numerical instability of adding floating point numbers,
  # the exact results of the summation dependent of the order of summation. For the purpose of generating reproducible numbers,
  # rounding the results to an arbitrary precision helps us achieve numerical stability.
  # This may be useful for unit-testing purposes
  if (!is.null(sigma.round)) results$sigma <- round(results$sigma, sigma.round)
  if (!is.null(beta.round)) results$beta <- lapply(results$beta, round, beta.round)
  
  results

}
