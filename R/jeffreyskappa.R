## Functions for Estimation (Sparse) Additive Topics (Kappa)
##########
#Contents: Kappa Estimation function, optimization of specific kappa, kappa likelihood/gradient.

# For the most complex case Kappa has 3 portions: 
#   topics, aspects and interactions (TAI).
# For Topics: we want to marginalize over the aspects.
# For Aspects: we want to marginalize over the topics.
# For Interaction: we just want the subset.
# In each case, c.k will be a vocab length vector.

jeffreysKappa <- function(beta.ss, kappa, settings) {
  verbose <- settings$verbose
  # Step 1: Precalulating some marginalizations
  if(length(beta.ss)==1) {
    KbyV <- beta.ss[[1]]
    KbyA <- matrix(rowSums(KbyV), ncol=1)
  } else {
    KbyV <- Reduce("+", beta.ss) #marginalize over aspects
    KbyA <- do.call(cbind, lapply(beta.ss, rowSums)) #marginalize over words
  }
  
  # Step 2: Loop over the Kappa Updates
  converged.outer <- FALSE
  outer.it <- 1
  
  while(!converged.outer) {
    compare <- abs(unlist(kappa$params)) <.001 #comparison for checking convergence at the end of the loop
    for(i in 1:length(kappa$params)) {
      kappa <- opt.kappak(i, kappa, beta.ss, 
                          settings$tau$mode, settings$tau$tol, settings$tau$maxit, settings$tau$prior,
                          KbyV, KbyA)
      if(verbose) cat(".")
    }
    if(verbose) cat("\n")  #end the update line
    
    #convergence checks for tolerance 
    current <- abs(unlist(kappa$params)) < .001
    sparseagree <- mean(compare==current)
    
    converged.outer <- sparseagree > settings$kappa$mstep$tol
    
    #convergence checks of iterations
    ifelse(outer.it==settings$kappa$mstep$maxit, converged.outer <- TRUE, outer.it <- outer.it + 1)
  }

  # Step 3: Construct Beta and Return 
  beta <- lapply(kappa$kappasum, function(x) exp(x - row.lse(x)))
  return(list(beta=beta, kappa=kappa))
}


opt.kappak <- function(i, kappa, beta.ss, taumode, tautol, taumaxit, taufixedprior, KbyV, KbyA) {
  a <- kappa$covar$a[i]
  k <- kappa$covar$k[i]
  type <- kappa$covar$type[i]
  
  kappa.init <- kappa$params[[i]]
  
  if(type==1) { #topics
    c.k <- KbyV[k,]
    bigC.k <- KbyA[k,]    
    kappa.other <- matrix(NA, nrow=ncol(KbyA), ncol=length(c.k))
    for(j in 1:ncol(KbyA)) {
      kappa$kappasum[[j]][k,] <- kappa$kappasum[[j]][k,] - kappa.init
      kappa.other[j,] <- exp(kappa$kappasum[[j]][k,,drop=FALSE])
    }
  }
  if(type==2) { #aspects
    c.k <- colSums(beta.ss[[a]])
    bigC.k <- KbyA[,a]
    kappa.other <- matrix(NA, nrow=nrow(KbyA), ncol=length(c.k))
    for(j in 1:nrow(KbyA)) {
      kappa$kappasum[[a]][j,] <- kappa$kappasum[[a]][j,] - kappa.init
      kappa.other[j,] <- exp(kappa$kappasum[[a]][j,,drop=FALSE])
    }
  }
  if(type==3) { #interactions
    c.k <- beta.ss[[a]][k,] 
    bigC.k <- KbyA[k,a]
    kappa$kappasum[[a]][k,] <- kappa$kappasum[[a]][k,] - kappa.init
    kappa.other <- exp(kappa$kappasum[[a]][k,,drop=FALSE])
  }
  
  converged <- FALSE
  its <- 1
  while(!converged) {
    gaussprec <- opt.tau(kappa.init, i, kappa, taumode, taufixedprior) 
    #recycle the vector for the c code
    if(length(gaussprec)!=length(c.k)) gaussprec <- rep(gaussprec, length(c.k))
    kappa.init <- optim(kappa.init, kappa.other=kappa.other, 
                        c.k=c.k, bigC.k=bigC.k,
                        gaussprec=gaussprec,
                        fn=kappa.obj, gr=kappa.gradient, 
                        method="L-BFGS-B", control=list(fnscale=-1, maxit=1000))$par
    #convergence checks
    converged <- mean(abs(kappa$params[[i]] - kappa.init)) <tautol 
    ifelse(its==taumaxit, converged <- TRUE, its <- its + 1)
    kappa$params[[i]] <- kappa.init
  }
  #add the bit back to kappa sum.
  if(type==1) { #topics, thus we need to add it to each aspect.
    for(j in 1:ncol(KbyA)) {
      kappa$kappasum[[j]][k,] <- kappa$kappasum[[j]][k,] + kappa$params[[i]]
    }  
  }
  if(type==2) { #aspects, thus we need to add to all topics in the aspect
    for(j in 1:nrow(kappa$kappasum[[a]])) {
      kappa$kappasum[[a]][j,] <- kappa$kappasum[[a]][j,] + kappa$params[[i]]
    }
  }
  if(type==3) { #interactions, only need to update in the one place!
    kappa$kappasum[[a]][k,] <- kappa$kappasum[[a]][k,] + kappa$params[[i]]
  }
  
  return(kappa)
}

###
# Optimization objectives for BFGS style Kappa
###
# in general we assume that kappa.other has the exp() 
# of the other kappa contributions precalculated.
# we have optimized these for speed with optim but further optimization
# would be possible in a setting where computations could be shared between
# the objective and gradient.
  
#the objective function is simply the sum of:
#(1) the kappa vector times the counts
#(2) the token counts with their appropriate denominator
#(3) the penalty
#it returns a scalar
kappa.obj <- function(kappa.param, kappa.other, c.k, bigC.k, gaussprec) {
  
  #(1) Kappa vector times counts
  p1 <- sum(c.k*kappa.param)
  
  #(2) Token Counts with Appropriate Denominator
  denom.kappas  <- sweep(kappa.other, MARGIN=2, exp(kappa.param), FUN="*")
  lseout <- log(rowSums(denom.kappas))
  p2 <- -sum(bigC.k*lseout)
  
  #(3) Sum over the penalty
  p3 <- -.5*sum(kappa.param^2*gaussprec)
  
  return((p1 + p2 + p3))
}

#the gradient has a similar interpretation:
#it is the difference of the observed counts and expected counts
#minus the penalty on kappa (divergence from zero, scaled by precision)
kappa.gradient <- function(kappa.param, kappa.other, c.k, bigC.k, gaussprec) {
  
  denom.kappas  <- sweep(kappa.other, MARGIN=2, exp(kappa.param), FUN="*")
  betaout <- denom.kappas/rowSums(denom.kappas)
  p2 <- -colSums(bigC.k*betaout)
    
  p3 <- -kappa.param*gaussprec
  
  return((c.k + p2 + p3))
}

