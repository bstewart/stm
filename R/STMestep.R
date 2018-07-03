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
                       verbose) {
  
  #quickly define useful constants
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  ctevery <- ifelse(N>100, floor(N/100), 1)
  if(!update.mu) mu.i <- as.numeric(mu)
  
  # 1) Initialize Sufficient Statistics 
  sigma.ss <- diag(0, nrow=(K-1))
  beta.ss <- vector(mode="list", length=A)
  for(i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
  
  # 2) Precalculate common components
  sigobj <- try(chol.default(sigma), silent=TRUE)
  if(class(sigobj)=="try-error") {
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
    sigmaentropy <- sum(log(diag(sigobj)))
    siginv <- chol2inv(sigobj)
  }
  # 3) Document Scheduling
  # For right now we are just doing everything in serial.
  # the challenge with multicore is efficient scheduling while
  # maintaining a small dimension for the sufficient statistics.
  for(i in 1:N) {
    #update components
    doc <- documents[[i]]
    words <- doc[1,]
    aspect <- beta.index[i]
    init <- lambda.old[i,]
    if(update.mu) mu.i <- mu[,i]
    beta.i <- beta[[aspect]][,words,drop=FALSE]
    
    #infer the document
    doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, 
                                  doc=doc, sigmaentropy=sigmaentropy)
    
    # update sufficient statistics 
    sigma.ss <- sigma.ss + doc.results$eta$nu
    beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
    bound[i] <- doc.results$bound
    lambda[[i]] <- c(doc.results$eta$lambda)
    if(verbose && i%%ctevery==0) cat(".")
  }
  if(verbose) cat("\n") #add a line break for the next message.
  
  #4) Combine and Return Sufficient Statistics
  lambda <- do.call(rbind, lambda)
  return(list(sigma=sigma.ss, beta=beta.ss, bound=bound, lambda=lambda))
}


#Let's start by assuming its one beta and we may have arbitrarily subset the number of docs.
estepmc <- function(documents, beta.index, update.mu, #null allows for intercept only model  
                    beta, lambda.old, mu, sigma, cores, 
                    verbose) {
  
  #quickly define useful constants
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  if(!update.mu) mu.i <- as.numeric(mu)
  
  # 1) Precalculate common components
  sigobj <- try(chol.default(sigma), silent=TRUE)
  if(class(sigobj)=="try-error") {
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
    sigmaentropy <- sum(log(diag(sigobj)))
    siginv <- chol2inv(sigobj)
  }

  groups <- base::split(seq_len(N),
                        sample(rep(seq_len(cores), length=N)))
  
  sysname <- Sys.info()["sysname"]
  if (sysname == "Windows") {
    cl <- parallel::makeCluster(cores)
    out <- parallel::clusterApply(cl, groups, fun=estep_parallel_block,
                                  documents=documents, beta.index=beta.index,
                                  lambda.old=lambda.old, mu=mu, beta=beta, siginv=siginv,
                                  sigmaentropy=sigmaentropy, update.mu=update.mu,
                                  K=K, N=N, V=V, A=A)
    parallel::stopCluster(cl)
  } else {
    out <- parallel::mclapply(groups, FUN=estep_parallel_block, mc.cores=cores,
                              documents=documents, beta.index=beta.index,
                              lambda.old=lambda.old, mu=mu, beta=beta, siginv=siginv,
                              sigmaentropy=sigmaentropy, update.mu=update.mu,
                              K=K, N=N, V=V, A=A)
  }
  
  beta.ss <- out[[1]]$beta.ss
  sigma.ss <- out[[1]]$sigma.ss
  
  for(i in 2:cores) {
    for(j in 1:A) {
      beta.ss[[j]] <- beta.ss[[j]] + out[[i]]$beta.ss[[j]]
    }
    sigma.ss <- sigma.ss + out[[i]]$sigma.ss
  }
  
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
  for(i in 1:cores) {
    bound[groups[[i]]] <- out[[i]]$bound
    lambda[groups[[i]]] <- out[[i]]$lambda
  }

  #4) Combine and Return Sufficient Statistics
  lambda <- do.call(rbind, lambda)
  return(list(sigma=sigma.ss, beta=beta.ss, bound=bound, lambda=lambda))
}

estep_parallel_block <- function(group, documents, beta.index, lambda.old,
                                 mu, beta, siginv, sigmaentropy, update.mu,
                                 K, N, V, A) {
  N <- length(group) #redefine N to group length
  
  # 1) Initialize Sufficient Statistics 
  sigma.ss <- diag(0, nrow=(K-1))
  beta.ss <- vector(mode="list", length=A)
  for(i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
  
  if(!update.mu) mu.i <- as.numeric(mu)
  
  # 2) Loop the documents in the group
  ct <- 0
  for(i in group) {
    ct <- ct + 1
    #update components
    doc <- documents[[i]]
    words <- doc[1,]
    aspect <- beta.index[i]
    init <- lambda.old[i,]
    if(update.mu) mu.i <- mu[,i]
    beta.i <- beta[[aspect]][,words,drop=FALSE]
    
    #infer the document
    doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, 
                                     doc=doc, sigmaentropy=sigmaentropy)
    
    # update sufficient statistics 
    sigma.ss <- sigma.ss + doc.results$eta$nu
    beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
    bound[ct] <- doc.results$bound
    lambda[[ct]] <- c(doc.results$eta$lambda)
  }
  return(list(sigma.ss=sigma.ss, beta.ss=beta.ss, bound=bound, lambda=lambda))
}