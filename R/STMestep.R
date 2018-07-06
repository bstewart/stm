#####################################
# SERIAL
#####################################
estepSerial <- function(N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv, verbose, cores=1) {
  
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
  if(verbose) cat("\n")
  
  lambda <- do.call(rbind, lambda)
  list(sigma=sigma.ss, beta=beta.ss, bound=bound, lambda=lambda)
}
#####################################

#####################################
# PARALLEL
#####################################
combineFn <- function(R, r) {
  R$sigma.ss <- R$sigma.ss + r$doc.results$eta$nu
  R$beta.ss[[r$aspect]][, r$words] <- R$beta.ss[[r$aspect]][, r$words] + r$doc.results$phis
  R$bound[r$i] <- r$doc.results$bound
  R$lambda[[r$i]] <- c(r$doc.results$eta$lambda)
  R
}

estepParallel <- function(N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv, verbose, cores=1) {
  
  # INITIALIZATION OF COMBINED RESULTS
  beta.ss <- vector(mode="list", length=A)
  for (i in 1:A) beta.ss[[i]] <- matrix(0, nrow=K, ncol=V)
  initt <- list(
    sigma.ss = diag(0, nrow=(K-1)),
    beta.ss = beta.ss,
    bound = vector(length=N),
    lambda = vector("list", length=N)
  )
  
  if (verbose) cat("Starting Parallel E-Step\n")
  
  res <- foreach (i = 1:N, .combine = combineFn, .multicombine = FALSE, .init = initt) %dopar% {
    doc <- documents[[i]]
    words <- doc[1,]
    aspect <- beta.index[i]
    init <- lambda.old[i,]
    if (update.mu) {
      mu.i <- mu[, i]
    } else {
      mu.i <- as.numeric(mu)
    }
    beta.i <- beta[[aspect]][, words, drop=FALSE]
    
    doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, doc=doc, sigmaentropy=sigmaentropy)
    list(i=i, doc.results=doc.results, aspect=aspect, words=words)
  }
  
  lambda <- do.call(rbind, res$lambda)
  list(sigma=res$sigma.ss, beta=res$beta.ss, bound=res$bound, lambda=lambda)
}
#####################################

estep <- function(documents, beta.index, update.mu, beta, lambda.old, mu, sigma, verbose, cores=1) {
  
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  if (!update.mu) mu.i <- as.numeric(mu)
  
  sigobj <- try(chol.default(sigma), silent=TRUE)
  if (class(sigobj)=="try-error") {
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
    sigmaentropy <- sum(log(diag(sigobj)))
    siginv <- chol2inv(sigobj)
  }
  
  f <- ifelse(cores>1, estepParallel, estepSerial)
  f(N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv, verbose, cores)
}
