#' Compute the gram matrix
#'
#' Take a Matrix object and compute a gram matrix
#'
#' Due to numerical error in floating points you can occasionally get
#' very small negative values.  Thus we check if the minimum is under 0 
#' and if true, assign elements less than 0 to zero exactly.  This is mostly
#' so we don't risk numerical instability later by introducing negative numbers of 
#' any sort.
#' @param mat a Matrix sparse Document by Term matrix
#' @export

gram <- function(mat) {
  nd <- rowSums(mat)
  mat <- mat[nd>=2,] #its undefined if we don't have docs of length 2
  nd <- nd[nd>=2]
  divisor <- nd*(nd-1)
  Htilde <- mat/sqrt(divisor)
  Hhat <- diag(colSums(mat/divisor))
  Q <- crossprod(Htilde) - Hhat
  #if(min(Q)<0) Q@x[Q@x < 0 ] <- 0
  return(as.matrix(Q))
}


#' Find Anchor Words
#'
#' Take a gram matrix Q and returns K anchors.
#' 
#' Does not use any randomness. Notes that anchors are such that you can request
#' more than you want but its not particularly easy to simply start in the middle
#' of a process.  One possibility for allowing this would be to return Qbar but
#' I'm not going to worry about that right now.
#'
#' @param Q the gram matrix
#' @param K the number of desired anchors
#' @param verbose if TRUE prints a dot to the screen after each anchor
#' @export
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

#' RecoverL2
#'
#' Recover the topic-word parameters from a set of anchor words using the RecoverL2
#' procedure of Arora et. al.
#' 
#' Using the exponentiated algorithm and an L2 loss identify the optimal convex 
#' combination of the anchor words which can reconstruct each additional word in the 
#' matrix.  Transform and return as a beta matrix.
#'
#' @param Q the gram matrix
#' @param anchor a vector of indices for rows of Q containing anchors
#' @param wprob the empirical word probabilities used to renorm the mixture weights. 
#' @param verbose if TRUE prints information as it progresses.
#' @param ... optional arguments that will be passed to the exponentiated gradient algorithm.
#' @return 
#' \item{A}{a matrix of dimension K by V.  This is acturally the transpose of A in Arora et al. and the matrix we call beta.}
#' \item{R}{a matrix of dimensions K by K that contains the topic covariances.}
#' \item{condprob}{a list of exponentiated gradient results.  useful for checking convergence.}
#' @export
recoverL2 <- function(Q, anchor, wprob, verbose=TRUE, ...) {
  Qbar <- Q/rowSums(Q)
  X <- Qbar[anchor,]
  XtX <- tcrossprod(X)
  
  #Word by Word Solve For the Convex Combination
  condprob <- vector(mode="list", length=nrow(Qbar))
  for(i in 1:nrow(Qbar)) {
    if(i %in% anchor) { 
      #if its an anchor we create a dummy entry that is 1 by definition
      vec <- rep(0, nrow(XtX))
      vec[match(i,anchor)] <- 1
      condprob[[i]] <- list(par=vec)
    } else {
      y <- Qbar[i,]
      condprob[[i]] <- expgrad(X,y,XtX, ...)
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
  weights <- lapply(condprob, function(x) x$par)
  weights <- do.call(rbind, weights)
  A <- weights*wprob
  A <- t(A)/colSums(A)
  
  #Recover The Topic-Topic Covariance Matrix
  Adag <- mpinv(A)  
  R <- t(Adag)%*%Q%*%Adag
  
  return(list(A=A, R=R, condprob=condprob))
}

#' Exponentiated Gradient with L2 Loss
#'
#' Find the optimal convex combination of features X which approximate the vector y under
#' and L2 loss.
#' 
#' An implementation of RecoverL2 based on David Mimno's code.  In this setting the
#' objective and gradient can both be kernalized leading to faster computations than
#' possible under the KL loss.  Specifically the computations are no longer dependent
#' on the size of the vocabulary making the operation essentially linear on the total
#' vocab size.
#' The matrix X'X can be passed as an argument in settings
#' like spectral algorithms where it is constant across multiple optimizations.  
#'
#' @param X the transposed feature matrix.
#' @param y the target vector
#' @param XtX optionally a precalculated crossproduct.
#' @param alpha an optional initialization of the parameter.  Otherwise it starts at 1/nrow(X).
#' @param tol convergence tolerance
#' @param max.its maximum iterations to run irrespective of tolerance
#' @return 
#' \item{par}{optimal weights}
#' \item{its}{number of iterations run}
#' \item{converged}{logical indicating if it converged}
#' \item{entropy}{entropy of the resulting weights}
#' \item{log.sse}{log of the sum of squared error}
#' @export
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