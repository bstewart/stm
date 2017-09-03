####
# Non-Exported Utility Functions
####

##
# Random Utilities

#function for collapsing a character vector to a comma separated list
#note that we eliminate things with zero length so that we can return
#fewer than n words and still have the lists look nice
commas <- function(text){  
  paste(text[nchar(text)>0], collapse=", ")
}

#from the R documentation for is.integer
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
} 

posint <- function(x) {
  all(is.wholenumber(x)) & all(x>0)
}


#' Draw from a Multivariate Normal
#' 
#' A basic function for doing multivariate normal simulations
#' via the cholesky decomposition of the covariance matrix. Function
#' is based on one by Peter Hoff.
#' 
#' This is a pretty standard multivariate normal generator. It could
#' almost certainly be faster if we ported it over to \pkg{RcppArmadillo}
#' but it isn't used a ton at the moment.
#' 
#' @param n number of draws
#' @param mu the K-dimensional mean
#' @param Sigma the K by K dimensional positive definite covariance matrix
#' @param chol.Sigma the cholesky decomposition of the Sigma matrix.
#' 
#' @keywords internal
rmvnorm<-function(n,mu,Sigma,chol.Sigma=chol(Sigma)) {
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol.Sigma) +c(mu))
}

#' Calculate FREX (FRequency and EXclusivity) Words
#' 
#' A primarily internal function for calculating FREX words.
#' We expect most users will use \code{\link{labelTopics}} instead.
#' 
#' FREX attempts to find words which are both frequent in and exclusive to a topic of interest.
#' Balancing these two traits is important as frequent words are often by themselves simply functional
#' words necessary to discuss any topic.  While completely exclusive words can be so rare as to not
#' be informative. This accords with a long-running trend in natural language processing which is best exemplified
#' by the Term frequency-Inverse document frequency metric.  
#' 
#' Our notion of FREX comes from a paper by Bischof and Airoldi (2012) which proposed a Hierarchical
#' Poisson Deconvolution model.  It relies on a known hierarchical structure in the documents and requires
#' a rather complicated estimation scheme.  We wanted a metric that would capture their core insight but
#' still be fast to compute.
#' 
#' Bischof and Airoldi consider as asummary for a word's contribution to a topic the harmonic mean of the
#' word's rank in terms of exclusivity and frequency.  The harmonic mean is attractive here because it 
#' does not allow a high rank along one of the dimensions to compensate for the lower rank in another. Thus
#' words with a high score must be high along both dimensions.
#' 
#' The formula is ' 
#'\deqn{FREX = \left(\frac{w}{F} + \frac{1-w}{E}\right)^{-1}}{FREX = ((w/F) + ((1-w)/E))^-1} 
#' where F is the frequency score given by the emperical CDF of the word in it's topic distribution.  Exclusivity
#' is calculated by column-normalizing the beta matrix (thus representing the conditional probability of seeing
#' the topic given the word).  Then the empirical CDF of the word is computed within the topic.  Thus words with
#' high values are those where most of the mass for that word is assigned to the given topic.
#' 
#' For rare words exclusivity will always be very high because there simply aren't many instances of the word.
#' If \code{wordcounts} are passed, the function will calculate a regularized form of this distribution using a
#' James-Stein type estimator described in \code{\link{js.estimate}}.
#' 
#' @param logbeta a K by V matrix containing the log probabilities of seeing word v conditional on topic k
#' @param w a value between 0 and 1 indicating the proportion of the weight assigned to frequency 
#' @param wordcounts a vector of word counts.  If provided, a James-Stein type shrinkage estimator is 
#' applied to stabilize the exclusivity probabilities. This helps with the concern that the rarest words
#' will always be completely exclusive.
#' @references 
#' Bischof and Airoldi (2012) "Summarizing topical content with word frequency and exclusivity"
#' In Proceedings of the International Conference on Machine Learning.
#' @seealso \code{\link{labelTopics}} \code{\link{js.estimate}}
#' @export
#' @keywords internal
calcfrex <- function(logbeta, w=.5, wordcounts=NULL) {
  excl <- t(t(logbeta) - col.lse(logbeta))
  if(!is.null(wordcounts)) {
    #if word counts provided calculate the shrinkage estimator
    excl <- safelog(sapply(1:ncol(excl), function(x) js.estimate(exp(excl[,x]), wordcounts[x])))
  } 
  freqscore <- apply(logbeta,1,data.table::frank)/ncol(logbeta)
  exclscore <- apply(excl,1,data.table::frank)/ncol(logbeta)
  frex <- 1/(w/freqscore + (1-w)/exclscore)
  apply(frex,2,order,decreasing=TRUE)
}

#' A James-Stein Estimator Shrinking to a Uniform Distribution
#' 
#' A primarily internal function used in \code{\link{calcfrex}}.
#' 
#' This calculates a James-Stein type shrinkage estimator for a discrete probability
#' distribution regularizing towards a uniform distribution. The amount of shrinkage
#' is a function of the variance of MLE and the L2 norm distance from the uniform.
#' 
#' This function is based off the ideas in Hausser and Strimmer (2009)
#' 
#' @param prob the MLE estimate of the discrete probability distribution
#' @param ct the count of words observed to estimate that distirbution
#' 
#' @references 
#' Hausser, Jean, and Korbinian Strimmer. "Entropy inference and the James-Stein estimator, 
#' with application to nonlinear gene association networks." Journal of Machine Learning Research 
#' 10.Jul (2009): 1469-1484.
#' @export
#' @keywords internal
js.estimate <- function(prob, ct) {
  if(ct<=1) {
    #basically if we only observe a count of 1
    #the variance goes to infinity and we get the uniform distribution.
    return(rep(1/length(prob), length(prob)))
  }
  # MLE of prob estimate
  mlvar <- prob*(1-prob)/(ct-1)
  unif <- rep(1/length(prob), length(prob)) 
  
  # Deviation from uniform
  deviation <- sum((prob-unif)^2)
  
  #take care of special case,if no difference it doesn't matter
  if(deviation==0) return(prob)
  
  lambda <- sum(mlvar)/deviation
  #if despite  our best efforts we ended up with an NaN number-just return the uniform distribution.
  if(is.nan(lambda)) return(unif)
  
  #truncate
  if(lambda>1) lambda <- 1
  if(lambda<0) lambda <- 0
  
  #Construct shrinkage estimator as convex combination of the two
  lambda*unif + (1 - lambda)*prob
}

#' Calculate Lift Words
#' 
#' A primarily internal function for calculating words according to the lift metric.
#' We expect most users will use \code{\link{labelTopics}} instead.
#' 
#' Lift is the calculated by dividing the topic-word distribution by the empirical
#' word count probability distribution.  In other words the Lift for word v in topic
#' k can be calculated as:
#' 
#' \deqn{Lift = \beta_{k,v}/(w_v/\sum_v w_v)}{Lift = \beta/wbar} 
#' 
#' We include this after seeing it used effectively in Matt Taddy's work including his
#' excellent \pkg{maptpx} package. Definitions are given in Taddy(2012).
#' 
#' @param logbeta a K by V matrix containing the log probabilities of seeing word v conditional on topic k
#' @param wordcounts a V length vector indicating the number of times each word appears in the corpus. 
#' @references 
#' Taddy, Matthew. 2012. "On Estimation and Selection for Topic Models." AISTATS JMLR W&CP 22
#' 
#' @seealso \code{\link{labelTopics}}
#' @export
#' @keywords internal
calclift <- function(logbeta, wordcounts) {
  emp.prob <- log(wordcounts) - log(sum(wordcounts))
  lift <- logbeta - rep(emp.prob, each=nrow(logbeta)) 
  apply(lift, 1, order, decreasing=TRUE)
}

#' Calculate Score Words
#' 
#' A primarily internal function for calculating words according to the score metric.
#' We expect most users will use \code{\link{labelTopics}} instead.
#' 
#' Score is a metric which we include because it is used effectively in the 
#' \pkg{lda} package by Jonathan Chang. It is calculated as:
#' \deqn{\beta_{v, k} (\log \beta_{w,k} - 1 / K \sum_{k'} \log \beta_{v,k'})}
#' 
#' @param logbeta a K by V matrix containing the log probabilities of seeing word v conditional on topic k
#' @references 
#' Jonathan Chang (2015). lda: Collapsed Gibbs Sampling Methods for Topic Models. R package version 1.4.2.
#' https://CRAN.R-project.org/package=lda 
#' @seealso \code{\link{labelTopics}} 
#' @export
#' @keywords internal
calcscore <- function(logbeta) { 
  ldascore <- exp(logbeta)*(logbeta - rep(colMeans(logbeta), each=nrow(logbeta)))
  apply(ldascore, 1, order, decreasing=TRUE)
} 
##
#Document convertors
##

#Our Format to Triplet Format
doc.to.ijv <- function(documents, fixzeroindex=TRUE) {
  #Turns our format into triplet format (can be zero indexed)
  indices <- unlist(lapply(documents, '[',1,)) #grab the first row
  if((0 %in% indices) & fixzeroindex) indices <- indices + 1 #if zero-indexed, fix it.
  counts <- lapply(documents, '[',2,)  #grab the second row but preserve the list structure for a moment
  VsubD <- unlist(lapply(counts,length)) #grab the number of unique words per document
  rowsums <- unlist(lapply(counts,sum)) #grab the number of tokens per documents
  docids <- rep(1:length(documents), times=VsubD) #add row numbers
  counts <- unlist(counts) #unlist the count structure
  #now we return using the convention for the simple triplet matrix,
  #plus the row sums which we use in DMR.
  return(list(i=as.integer(docids), j=as.integer(indices), v=as.integer(counts), rowsums=as.integer(rowsums)))
}

#Triplet Format to our Document Format
ijv.to.doc <- function(i,j,v) {
  index <- split(j,i)
  index <- lapply(index,as.integer)
  count <- split(v,i)
  count <- lapply(count,as.integer)
  mapply(rbind,index,count)
}

##
#A series of fast softmax functions mostly wrappers around matrixStats package functions
##
logsoftmax <- function(x) {
  x - lse(x)
}

lse <- function(x) {
 matrixStats::logSumExp(x)
}

row.lse <- function(mat) {
  matrixStats::rowLogSumExps(mat)
}
col.lse <- function(mat) {
  matrixStats::colLogSumExps(mat)
}

softmax <- function(x) {
  exp(x - lse(x))
}

safelog <- function(x, min=-1000) {
  out <- log(x)
  out[which(out< min)] <- min
  out
}


# Note: I started documenting this but am not exporting because I would need to appropriately
# generalize.  It is currently only set up to do the mgaussian one I think- but I haven't looked
# at it in a while.

#' Unpack a \pkg{glmnet} object
#' 
#' A function to quickly unpack a \pkg{glmnet} model object and calculate an
#' optimal model from the regularization path.
#' 
#' This is a small utility we wrote to deal with the slow methods dispatch for S4
#' classes.  The more straightforward option is the \code{coef()} method for \pkg{glmnet}
#' objects but when trying to make thousands of calls a second, that can be very slow
#' 
#' @param mod the glmnet model
#' @param ic.k the information criterion value.  AIC is \code{ic.k=2} and BIC would be \code{ic.k=log n}
#' 
#' @return 
#' A list
#' \item{coef}{a matrix of coefficients}
#' \item{intercept}{the intercepts}
#' @keywords internal
unpack.glmnet <- function(mod, ic.k) {
  dev <- (1-mod$dev.ratio)*mod$nulldev
  df <- colSums(mod$dfmat)
  ic <- dev + ic.k*df
  lambda <- which.min(ic)
  
  #methods dispatch here is crazy, so define out own function
  subM <- function(x, p) {
    ind <- (x@p[p]+1):x@p[p+1]
    rn <- x@i[ind]+1
    y <- x@x[ind]
    out <- rep(0, length=nrow(x))
    out[rn] <- y
    out
  }
  coef <- lapply(mod$beta, subM, lambda) #grab non-zero coefs
  coef <- do.call(cbind,coef)
  intercept <- mod$a0[,lambda]
  return(list(coef=coef, intercept=intercept))
}
