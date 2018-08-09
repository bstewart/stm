#####################################
# SERIAL
#####################################
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
  if(verbose) cat("\n")
  
  lambda <- do.call(rbind, lambda)
  list(sigma=sigma.ss, beta=beta.ss, bound=bound, lambda=lambda)
}
#####################################

#####################################
# PARALLEL
#####################################
combineFn <- function(R, r) {
  R$sigma.ss <- R$sigma.ss + r$sigma.ss
  for (i in length(R$beta.ss)) {
    R$beta.ss[[i]] =  R$beta.ss[[i]] + r$beta.ss[[i]]
  }
  R$bound[r$doc.ids] <- r$bound[r$doc.ids]
  R$lambda[r$doc.ids] <- r$lambda[r$doc.ids]
  R
}

#' @importFrom foreach foreach %dopar%
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
  
  # doc.id.groups <- base::split(seq_len(N), sample(rep(seq_len(cores), length=N)))
  doc.id.groups <- base::split(seq_len(N), rep(seq_len(cores), length=N))
  
  res <- foreach (doc.ids = doc.id.groups, .combine = combineFn, .multicombine = FALSE, .init = initt) %dopar% {
    estepParallelBlock(doc.ids, N, K, A, V, documents[doc.ids], beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv)
  }
  
  lambda <- do.call(rbind, res$lambda)
  list(sigma=res$sigma.ss, beta=res$beta.ss, bound=res$bound, lambda=lambda)
}

estepParallelBlock <- function(doc.ids, N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv) {
  
  sigma.ss <- diag(0, nrow=K-1)
  beta.ss <- vector(mode='list', length=A)
  for(i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow=K, ncol=V)
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)

  if(!update.mu) mu.i <- as.numeric(mu)
  
  cnt <- 1
  for (i in doc.ids) {
    doc = documents[[cnt]]
    words <- doc[1,]
    aspect <- beta.index[i]
    init <- lambda.old[i,]
    if (update.mu) mu.i <- mu[, i]
    beta.i <- beta[[aspect]][, words, drop=FALSE]
    
    doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, doc=doc, sigmaentropy=sigmaentropy)
    
    sigma.ss <- sigma.ss + doc.results$eta$nu
    beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
    bound[i] <- doc.results$bound
    lambda[[i]] <- c(doc.results$eta$lambda)
    
    cnt <- cnt + 1
    
  }
  list(doc.ids=doc.ids, sigma.ss=sigma.ss, beta.ss=beta.ss, bound=bound, lambda=lambda)
}

#####################################

estep <- function(documents, beta.index, update.mu, beta, lambda.old, mu, sigma, verbose, cores=1, round.estep.sigma.digits=NULL, round.estep.beta.digits=NULL) {
  
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
  
  if (cores>1) {
    results <- estepParallel(N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv, verbose, cores)
  } else {
    results <- estepSerial(N, K, A, V, documents, beta.index, lambda.old, mu, update.mu, beta, sigmaentropy, siginv, verbose)
  }
  
  if (!is.null(round.estep.sigma.digits)) results$sigma <- round(results$sigma, round.estep.sigma.digits)
  if (!is.null(round.estep.beta.digits)) results$beta <- lapply(results$beta, round, round.estep.beta.digits)
  
  results

}
