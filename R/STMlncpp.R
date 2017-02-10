
#The internal unexported function (external wrapper is below).
#incidentally for those monkeying around down here in the code
#the wrapper provides some useful documentation.  Essentially all
#it does is subset beta for the user and make the argument ordering
#and naming matching the rest of the package.

#12/31/2016 added the hpbcpp argument so I can opt not to call it
#inside fitNewDocuments
logisticnormalcpp <- function(eta, mu, siginv, beta, doc, sigmaentropy, 
                              method="BFGS", control=list(maxit=500),
                              hpbcpp=TRUE) {
  doc.ct <- doc[2,]
  Ndoc <- sum(doc.ct)
  #even at K=100, BFGS is faster than L-BFGS
  optim.out <- optim(par=eta, fn=lhoodcpp, gr=gradcpp,
                     method=method, control=control,
                     doc_ct=doc.ct, mu=mu,
                     siginv=siginv, beta=beta)
  
  if(!hpbcpp) return(list(eta=list(lambda=optim.out$par)))
  
  #Solve for Hessian/Phi/Bound returning the result
  hpbcpp(optim.out$par, doc_ct=doc.ct, mu=mu,
         siginv=siginv, beta=beta,
         sigmaentropy=sigmaentropy)
}


#the external exported version
####

#' Optimize Document
#' 
#' A primarily internal use function for optimizing the document-level
#' parameters of the variational distribution.  
#' Included here for advanced users who want to design new
#' post-processing features. This help file assumes knowledge of our 
#' notation which follows the mathematical notation used in our vignette
#' and other papers.
#' 
#' This function is a small wrapper around the internal function used
#' to complete the E-step for each document.
#' 
#' Regarding the arguments \code{sigma}, \code{sigmainv} and \code{sigmaentropy}.  In
#' the internal version of the code we calculate \code{sigmainv} and \code{sigmaentropy}
#' once each E-step because it is shared by all documents.  If you supply the original
#' value to \code{sigma} it will calculate these for you.  If you are going to be using
#' this to run a bunch of documents and speed is a concern, peek at the underlying code
#' and do the calculation youself once and then just pass the result to the function so
#' it isn't repeated with every observation.  
#' 
#' @param document a single matrix containing the document in the \code{\link{stm}} format
#' @param eta a vector of length K-1 containing the initial starting value for eta
#' @param mu a vector of length K-1 containing the prevalence prior
#' @param beta a matrix containing the complete topic-word distribution for the document.
#' If using a content covariate model it is presumed that you have already passed the correct content
#' covariate level's beta.
#' @param sigma a K-1 by K-1 matrix containing the covariance matrix of the MVN prior. If you supply this
#' you do not need to supply \code{sigmainv} or \code{sigmaentropy}.  See below.
#' @param sigmainv a K-1 by K-1 matrix containing the precision matrix of the MVN prior.  If you supplied
#' \code{sigma} you do not need to supply this. See below.
#' @param sigmaentropy the entropy term calculated from sigma.  If you supplied \code{sigma} you do not
#' need to supply this.  See below.
#' @param method the method passed to \code{\link{optim}}.  Uses "BFGS" by default.
#' @param control the control argument passed to \code{\link{optim}}.  Sets the maximum number of observations
#' to 500 but can be used to set other aspects of the optimization per the instructions in \code{\link{optim}}
#' @param posterior should the full posterior be returned?  If TRUE (as it is by default) returns the full 
#' variational posterior.  Otherwise just returns the point estimate.
#' @return a list 
#' 
#' \item{phis}{A K by V* matrix containing the variational distribution for each token (where V* is the number of 
#' unique words in the given document.  They are in the order of appearence in the document. For words repeated
#' more than once the sum of the column is the number of times that token appeared.}
#' \item{lambda}{A (K-1) by 1 matrix containing the mean of the variational distribution for eta.  This is 
#' actually just called eta in the output of \code{\link{stm}} as it is also the point estimate.}
#' \item{nu}{A (K-1) by (K-1) matrix containg the covariance matrix of the variational distribution for eta.
#' This is also the inverse Hessian matrix.}
#' \item{bound}{The value of the document-level contribution to the global approximate evidence lower bound.}
#' 
#' @seealso \code{\link{thetaPosterior}} 
#' @examples
#' # fitting to a nonsense word distribution
#' V <- length(poliblog5k.voc)
#' K <- 50
#' beta <- matrix(rgamma(V*K,shape = .1), nrow=K, ncol=V)
#' beta <- beta/rowSums(beta)
#' doc <- poliblog5k.docs[[1]]
#' mu <- rep(0, K-1)
#' sigma <- diag(1000, nrow=K-1)
#' optimizeDocument(doc, eta=rep(0, K-1), mu=mu, beta=beta, sigma=sigma)
#' @export
optimizeDocument <- function(document, eta, mu, beta, sigma=NULL, 
                             sigmainv=NULL, sigmaentropy=NULL,
                             method="BFGS", control=list(maxit=500),
                             posterior=TRUE) {
  if(is.null(sigma) & (is.null(sigmainv) | is.null(sigmaentropy))) {
    stop("You must specify either sigma OR sigma inverse and sigma entropy")
  }
  if(!is.matrix(document)) {
    stop("document must be a single matrix. See help file.")
  }
  
  if(!is.null(sigma)) {
    sigobj <- try(chol.default(sigma), silent=TRUE)
    if(class(sigobj)=="try-error") {
      sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
      siginv <- solve(sigma)
    } else {
      sigmaentropy <- sum(log(diag(sigobj)))
      siginv <- chol2inv(sigobj)
    }
  } else {
    siginv <- sigmainv #just renaming otherwise
  }
  if(any(beta < 0)) stop("Some entries of beta are negative.  Are you sure you
                         didn't pass the logged version of beta?")
  
  beta <- beta[,document[1,]]
  out <- logisticnormalcpp(eta, mu, siginv, beta, document, sigmaentropy,
                           method=method, control=control, hpbcpp=posterior)
  toReturn <- list(lambda=out$eta$lambda, phi=out$phi, nu=out$eta$nu, bound=out$bound)
  return(toReturn)
}


 