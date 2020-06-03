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
                       summation, randomize, method,
                       verbose) {
  #quickly define useful constants
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  ctevery <- ifelse(N>100, floor(N/100), 1)
  if(!update.mu) mu.i <- as.numeric(mu)
  # 1) Initialize Sufficient Statistics 
  if(summation$neum_R) {
    sigma.ss <- n_mat_sum(diag(0, nrow=(K-1)))
  } else if(summation$neum_cpp) {
    sigma.ss <- matrix(0, nrow=2*(K-1), ncol=K-1)
  } else {
    sigma.ss <- diag(0, nrow=(K-1))
  }
  if(summation$neum_R) {
    beta.ss <- vector(mode="list", length=A)
    for(i in 1:A) {
      beta.ss[[i]] <- n_mat_sum(matrix(0, nrow=K,ncol=V))
    }
  } else if(summation$neum_cpp) {
    beta.ss <- vector(mode="list", length=A)
    for(i in 1:A) {
      beta.ss[[i]] <- matrix(0, nrow=2*K,ncol=V)
    }
  }
  else {
    beta.ss <- vector(mode="list", length=A)
    for(i in 1:A) {
      beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
    }
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
  if(randomize) {
    vec <- sample(1:N, N) 
  } else {
    vec <- 1:N
  }

  for(i in vec) {
    #update components
    doc <- documents[[i]]
    words <- doc[1,]
    aspect <- beta.index[i]
    init <- lambda.old[i,]
    if(update.mu) mu.i <- mu[,i]
    beta.i <- beta[[aspect]][,words,drop=FALSE]
    
    #infer the document
    doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, doc=doc, sigmaentropy=sigmaentropy, method=method)
    
    # update sufficient statistics 
    if(summation$neum_R) {
      sigma.ss <- n_mat_sum(sigma.ss[[1]], sigma.ss[[2]], doc.results$eta$nu)
    } else if(summation$neum_cpp) {
      n_sigma_sumcpp(sigma.ss, doc.results$eta$nu)
    } else {
      sigma.ss <- sigma.ss + doc.results$eta$nu
    }
    if(summation$neum_R) {
      #more efficient than this would be to stack all the C's underneath
      #betas
      o_beta <- n_mat_sum(beta.ss[[aspect]][[1]][,words], 
                          beta.ss[[aspect]][[2]][,words], 
                          doc.results$phis)
      beta.ss[[aspect]][[1]][,words] <- o_beta[[1]]
      beta.ss[[aspect]][[2]][,words] <- o_beta[[2]]
    } else if(summation$neum_cpp) {
      n_beta_sumcpp(beta.ss[[aspect]], words-1, doc.results$phis)
    } else {
      beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
    }
    bound[i] <- doc.results$bound
    lambda[[i]] <- c(doc.results$eta$lambda)
    if(verbose && i%%ctevery==0) cat(".")
  }
  if(verbose) cat("\n") #add a line break for the next message.
  
  #4) Combine and Return Sufficient Statistics
  lambda <- do.call(rbind, lambda)
  if(summation$neum_cpp) {
    sigma.ss <- matrix((t(bdiag(rep(list(matrix(c(1,1,0,1), nrow=2, ncol=2)), as.integer(nrow(sigma.ss)/2)))) %*% sigma.ss)[-seq(0, nrow(sigma.ss), 2),], nrow=as.integer(nrow(sigma.ss)/2), ncol=ncol(sigma.ss))
    beta.ss <- lapply(beta.ss, function(x) matrix((t(bdiag(rep(list(matrix(c(1,1,0,1), nrow=2, ncol=2)), as.integer(nrow(x)/2)))) %*% x)[-seq(0, nrow(x), 2),], nrow=as.integer(nrow(x)/2), ncol=ncol(x)))
  }
  return(list(sigma=sigma.ss, beta=beta.ss, bound=bound, lambda=lambda, vec=vec))
}
