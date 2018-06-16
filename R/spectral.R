# Compute the gram matrix
#
# Take a Matrix object and compute a gram matrix
#
# Due to numerical error in floating points you can occasionally get
# very small negative values.  Thus we check if the minimum is under 0 
# and if true, assign elements less than 0 to zero exactly.  This is mostly
# so we don't risk numerical instability later by introducing negative numbers of 
# any sort.
# @param mat A Matrix sparse Document by Term matrix
gram <- function(mat) {
  nd <- Matrix::rowSums(mat)
  mat <- mat[nd>=2,] #its undefined if we don't have docs of length 2
  nd <- nd[nd>=2]
  divisor <- nd*(nd-1)
  #clearer code is removed below for memory efficiency
  #Htilde <- mat/sqrt(divisor)
  #Hhat <- diag(colSums(mat/divisor))
  #Q <- crossprod(Htilde) - Hhat
  Q <- Matrix::crossprod(mat/sqrt(divisor)) - Matrix::diag(Matrix::colSums(mat/divisor))
  #if(min(Q)<0) Q@x[Q@x < 0 ] <- 0
  return(as.matrix(Q))
}

gram.rp <- function(mat, s=.05, p=3000, d.group.size=2000, verbose=TRUE) {
  # s is the sparsity level
  # p is the number of projections
  # d.group.size is the size of the groups of documents analyzed.
  
  #first, figure out how many items are necessary to produce the correct sparsity level
  n.items <- ceiling(s*p)
  #construct a projection matrix out of a triplet representation
  triplet <- vector(mode="list", length=ncol(mat))
  for(i in 1:ncol(mat)) {
    index <- sample(size=n.items, x=1:p, replace=FALSE)
    values <- sample(size=length(index), x=c(-1, 1), replace=TRUE)
    triplet[[i]] <- cbind(i, index, values)
  }
  triplet <- do.call(rbind,triplet)
  proj <- sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3])
  rm(triplet)
  
  #now we do the projection.
  # the basic strategy is to do the computation in blocks of documents
  # I tried various methods of doing one document at a time but they 
  # were all too slow.
  D <- nrow(mat)
  groups <- ceiling(D/d.group.size)
  Qnorm <- rep(0, ncol(mat))
  #iterate over blocks adding to total
  if(verbose) cat("\t")
  for(i in 1:groups) {
    smat <- mat[((i-1)*d.group.size + 1):min(i*d.group.size, D),]
    #not positive about next two lines.
    rsums <- rowSums(smat)
    divisor <- rsums*(rsums-1)
    #divide through
    smatd <- smat/divisor
    #update the piece that we will ultimately normalize with
    Qnorm <- Qnorm + colSums(smatd*(rsums - 1 - smatd))
    #calculate for the bit below
    Htilde <- colSums(smatd)*proj
    rm(smatd) #clear out a copy we don't need
    smat <- smat/sqrt(divisor)
    
    if(i>1) {
      Q <- Q + t(smat)%*%(smat%*%proj) - Htilde
    } else {
      #if its the first 
      Q <- t(smat)%*%(smat%*%proj) - Htilde
    }
    if(verbose) cat(".")
  }
  if(verbose) cat("\n")
  return(as.matrix(Q/Qnorm))
}

# Find Anchor Words
#
# Take a gram matrix Q and returns K anchors.
# 
# Does not use any randomness. Notes that anchors are such that you can request
# more than you want but its not particularly easy to simply start in the middle
# of a process.  One possibility for allowing this would be to return Qbar but
# I'm not going to worry about that right now.
#
# @param Q The gram matrix
# @param K The number of desired anchors
# @param verbose If TRUE prints a dot to the screen after each anchor
fastAnchor <- function(Qbar, K, verbose=TRUE) {
  basis <- c()
  rowSquaredSums <- rowSums(Qbar^2) #StabilizedGS
  
  for(i in 1:K) {
    basis[i] <- which.max(rowSquaredSums) #83-94
    
    maxval <- rowSquaredSums[basis[i]]
    normalizer <- 1/sqrt(maxval) #100
    
    #101-103
    Qbar[basis[i],] <- Qbar[basis[i],]*normalizer 
    
    #For each row
    innerproducts <- Qbar%*%Qbar[basis[i],] #109-113
    
    #Each row gets multiplied out through the Basis
    project <- as.numeric(innerproducts)%o%Qbar[basis[i],] #118
    
    #Now we want to subtract off the projection but
    #first we should zero out the components for the basis
    #vectors which we weren't intended to calculate
    project[basis,] <- 0 #106 (if statement excluding basis vectors)
    Qbar <- Qbar - project #119
    rowSquaredSums <- rowSums(Qbar^2)
    rowSquaredSums[basis] <- 0 #here we cancel out the components we weren't calculating.
    if(verbose) cat(".")
  }
  return(basis)
}

# RecoverL2
#
# Recover the topic-word parameters from a set of anchor words using the RecoverL2
# procedure of Arora et. al.
# 
# Using the exponentiated algorithm and an L2 loss identify the optimal convex 
# combination of the anchor words which can reconstruct each additional word in the 
# matrix.  Transform and return as a beta matrix.
#
# @param Qbar The row-normalized gram matrix
# @param anchor A vector of indices for rows of Q containing anchors
# @param wprob The empirical word probabilities used to renorm the mixture weights. 
# @param verbose If TRUE prints information as it progresses.
# @param ... Optional arguments that will be passed to the exponentiated gradient algorithm.
# @return 
# \item{A}{A matrix of dimension K by V.  This is acturally the transpose of A in Arora et al. and the matrix we call beta.}
recoverL2 <- function(Qbar, anchor, wprob, verbose=TRUE, recoverEG=TRUE, ...) {
  #NB: I've edited the script to remove some of the calculations by commenting them
  #out.  This allows us to store only one copy of Q which is more memory efficient.
  #documentation for other pieces is below.

  #Qbar <- Q/rowSums(Q)
  X <- Qbar[anchor,]
  XtX <- tcrossprod(X)
  
  #In a minute we will do quadratic programming
  #these jointly define the conditions.  First column
  #is a sum to 1 constraint.  Remainder are each parameter
  #greater than 0.
  Amat <- cbind(1,diag(1,nrow=nrow(X)))
  bvec <- c(1,rep(0,nrow(X)))
  
  #Word by Word Solve For the Convex Combination
  condprob <- vector(mode="list", length=nrow(Qbar))
  for(i in 1:nrow(Qbar)) {
    if(i %in% anchor) { 
      #if its an anchor we create a dummy entry that is 1 by definition
      vec <- rep(0, nrow(XtX))
      vec[match(i,anchor)] <- 1
      condprob[[i]] <- vec
    } else {
      y <- Qbar[i,]
      
      if(recoverEG) {
        solution <- expgrad(X,y,XtX, ...)$par
      } else {
        #meq=1 means the sum is treated as an exact equality constraint
        #and the remainder are >=
        solution <- quadprog::solve.QP(Dmat=XtX, dvec=X%*%y, 
                                       Amat=Amat, bvec=bvec, meq=1)$solution  
      }
      
      if(any(solution <= 0)) {
        #we can get exact 0's or even slightly negative numbers from quadprog
        #replace with machine double epsilon
        solution[solution<=0] <- .Machine$double.eps
      } 
      condprob[[i]] <- solution
    }
    if(verbose) {
      #if(i%%1 == 0) cat(".")
      #if(i%%20 == 0) cat("\n")
      #if(i%%100 == 0) cat(sprintf("Recovered %i of %i words. \n", i, nrow(Qbar)))
      if(i%%100==0) cat(".")
    }
  }
  if(verbose) cat("\n")
  #Recover Beta (A in this notation)
  #  Now we have p(z|w) but we want the inverse
  weights <- do.call(rbind, condprob)
  A <- weights*wprob
  A <- t(A)/colSums(A)
  
  #Recover The Topic-Topic Covariance Matrix
  #Adag <- mpinv(A)  
  #R <- t(Adag)%*%Q%*%Adag
  return(list(A=A))
  #return(list(A=A, R=R, condprob=condprob))
}

# Exponentiated Gradient with L2 Loss
#
# Find the optimal convex combination of features X which approximate the vector y under
# and L2 loss.
# 
# An implementation of RecoverL2 based on David Mimno's code.  In this setting the
# objective and gradient can both be kernalized leading to faster computations than
# possible under the KL loss.  Specifically the computations are no longer dependent
# on the size of the vocabulary making the operation essentially linear on the total
# vocab size.
# The matrix X'X can be passed as an argument in settings
# like spectral algorithms where it is constant across multiple optimizations.  
#
# @param X The transposed feature matrix.
# @param y The target vector
# @param XtX Optionally a precalculated crossproduct.
# @param alpha An optional initialization of the parameter.  Otherwise it starts at 1/nrow(X).
# @param tol Convergence tolerance
# @param max.its Maximum iterations to run irrespective of tolerance
# @return 
# \item{par}{Optimal weights}
# \item{its}{Number of iterations run}
# \item{converged}{Logical indicating if it converged}
# \item{entropy}{Entropy of the resulting weights}
# \item{log.sse}{Log of the sum of squared error}
expgrad <- function(X, y, XtX=NULL, alpha=NULL, tol=1e-7, max.its=500) {
  if(is.null(alpha)) alpha <- 1/nrow(X) 
  alpha <- matrix(alpha, nrow=1, ncol=nrow(X))
  if(is.null(XtX)) XtX <- tcrossprod(X)
  
  ytX <- y%*%t(X)
  converged <- FALSE
  eta <- 50
  sse.old <- Inf
  its <- 1
  while(!converged) {
    #find the gradient (y'X - alpha'X'X)
    grad <- (ytX - alpha%*%XtX) #101-105
    sse <- sum(grad^2) #106, sumSquaredError
    grad <- 2*eta*grad
    maxderiv <- max(grad)    
    
    #update parameter 
    alpha <- alpha*exp(grad-maxderiv)
    #project parameter back to space
    alpha <- alpha/sum(alpha)
    
    converged <- abs(sqrt(sse.old)-sqrt(sse)) < tol
    if(its==max.its) break
    sse.old <- sse
    its <- its + 1
  } 
  entropy <- -1*sum(alpha*log(alpha))
  return(list(par=as.numeric(alpha), its=its, converged=converged,
              entropy=entropy, log.sse=log(sse)))
}

#calculate the Moore Penrose Pseudo-Inverse
# adapted from the mASS function ginv
mpinv <- function (X) {

  if (!is.matrix(X)) X <- as.matrix(X)
  
  tol <- sqrt(.Machine$double.eps)
  Xsvd <- svd(X)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                                               t(Xsvd$u[, Positive, drop = FALSE]))
}

tsneAnchor <- function(Qbar, verbose=TRUE, init.dims=50, perplexity=30) {
  if(!(requireNamespace("Rtsne",quietly=TRUE) & requireNamespace("geometry", quietly=TRUE & requireNamespace("rsvd", quietly=TRUE)))){
    stop("Please install the Rtsne, rsvd and geometry packages to use this setting.")
  } 
  if(verbose) cat("\t Initializing tSNE with PCA...\n \t")
  Xpca <- rsvd::rpca(Qbar, min(init.dims,ncol(Qbar)), center=TRUE, scale=FALSE, retx=TRUE)$x[,1: min(init.dims,ncol(Qbar))]
  #project to 3-D
  if(verbose) cat("\t Using tSNE to project to a low-dimensional space...\n \t")
  proj <- try(Rtsne::Rtsne(Xpca, pca=FALSE, dims=3,
                           initial_dims=init.dims,
                           perplexity=perplexity) , silent=TRUE)
  if(class(proj)=="try-error") {
    if(attr(proj, "condition")$message=="Perplexity is too large.") {
      rate <- 1
      while(class(proj)=="try-error" && attr(proj, "condition")$message=="Perplexity is too large.") {
        rate <- rate + 1
        if(verbose) cat(sprintf("\t tSNE failed because perplexity is too large. Using perplexity=%f \n \t", perplexity/rate))
        proj <- try(Rtsne::Rtsne(Xpca, pca=FALSE, dims=3,
                                 initial_dims=init.dims,
                                 perplexity=perplexity/rate) , silent=TRUE)
      }
    } else {
      #if this failed it is probably duplicates which Rtsne cannot handle
      dup <- duplicated(Xpca)
      if(!any(dup)) stop("an unknown error has occured in Rtsne")
      
      dup <- which(dup)
      for(r in dup) {
        row <- Qbar[r,]
        row[row!=0] <- runif(sum(row!=0),0,1e-5) # add a bit of noise to non-zero duplicates
        row <- row/sum(row) #renormalize
        Qbar[r,] <- row
      }
      #we now have to reproject because we have messed with Qbar
      Xpca <- rsvd::rpca(Qbar, min(init.dims,ncol(Qbar)), center=TRUE, scale=FALSE, retx=TRUE)$x[,1: min(init.dims,ncol(Qbar))]
      
      #and now do it again
      proj <- Rtsne::Rtsne(Xpca, pca=FALSE, dims=3,
                           initial_dims=init.dims,
                           perplexity=perplexity)
    }
  }
  if(verbose) cat("\t Calculating exact convex hull...\n \t")
  hull <- geometry::convhulln(proj$Y)
  anchor <- sort(unique(c(hull)))
  return(anchor)
}
