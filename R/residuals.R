#' Residual dispersion test for topic number
#' 
#' Computes the multinomial dispersion of the STM residuals as in Taddy (2012)
#' 
#' This function implements the residual-based diagnostic method of Taddy
#' (2012).  The basic idea is that when the model is correctly specified the
#' multinomial likelihood implies a dispersion of the residuals:
#' \eqn{\sigma^2=1}.  If we calculate the sample dispersion and the value is
#' greater than one, this implies that the number of topics is set too low,
#' because the latent topics are not able to account for the overdispersion. In
#' practice this can be a very demanding criterion, especially if the documents
#' are long.  However, when coupled with other tools it can provide a valuable
#' perspective on model fit. The function is based on the Taddy 2012 paper as well as code
#' found in maptpx package.
#' 
#' Further details are available in the referenced paper, but broadly speaking
#' the dispersion is derived from the mean of the squared adjusted residuals.
#' We get the sample dispersion by dividing by the degrees of freedom
#' parameter.  In estimating the degrees of freedom, we follow Taddy (2012) in
#' approximating the parameter \eqn{\hat{N}} by the number of expected counts
#' exceeding a tolerance parameter.  The default value of 1/100 given in the
#' Taddy paper can be changed by setting the \code{tol} argument.
#' 
#' The function returns the estimated sample dispersion (which equals 1 under
#' the data generating process) and the p-value of a chi-squared test where the
#' null hypothesis is that \eqn{\sigma^2=1} vs the alternative \eqn{\sigma^2
#' >1}. As Taddy notes and we echo, rejection of the null 'provides a very
#' rough measure for evidence in favor of a larger number of topics.'
#' 
#' @param stmobj An \code{STM} model object for which to compute residuals.
#' @param documents The documents corresponding to \code{stmobj} as in
#' \code{\link{stm}}.
#' @param tol The tolerance parameter for calculating the degrees of freedom.
#' Defaults to 1/100 as in Taddy(2012)
#' @references Taddy, M. 'On Estimation and Selection for Topic Models'.
#' AISTATS 2012, JMLR W&CP 22
#' @examples
#' 
#' #An example using the Gadarian data.  From Raw text to fitted model.
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' meta<-temp$meta
#' vocab<-temp$vocab
#' docs<-temp$documents
#' out <- prepDocuments(docs, vocab, meta)
#' docs<-out$documents
#' vocab<-out$vocab
#' meta <-out$meta
#' set.seed(02138)
#' #maximum EM iterations set very low so example will run quickly.  
#' #Run your models to convergence!
#' mod.out <- stm(docs, vocab, 3, prevalence=~treatment + s(pid_rep), data=meta,
#'                max.em.its=5)
#' checkResiduals(mod.out, docs)
#' @export
checkResiduals <- function(stmobj, documents, tol=.01) {
  
  # Convert the corpus to the internal STM format
  args <- asSTMCorpus(documents)
  documents <- args$documents

  beta <- lapply(stmobj$beta$logbeta, exp)
  theta <- stmobj$theta
  index <- stmobj$settings$covariates$betaindex
  K <- stmobj$settings$dim$K
  n <- length(documents)
  phat <- stmobj$settings$dim$V
  d <- n*(K-1) + K*( phat-1 )
  
  doc.resid <- function(doc, theta, beta) {
    q <- theta%*%beta 
    m <- sum(doc[2,])
    Nhat <- sum(q*m > tol)
    x <- rep(0, ncol(beta))
    x[doc[1,]] <- doc[2,]
    out <- sum((x^2 - 2*x*q*m)/(m*q*(1-q))) + sum(m*q/(1-q))
    return(list(out=out, Nhat=Nhat))
  }
  result <- vector(mode="list", length=length(documents))
  Nhat <- 0
  for(i in 1:length(result)) {
    resid <- doc.resid(documents[[i]], theta[i,], beta[[index[i]]])
    result[[i]] <- resid$out
    Nhat <- Nhat + resid$Nhat
  }
  D <- sum(unlist(result))
  df <- Nhat - phat  - d
  sig2 <- D/df
  
  rho <- suppressWarnings(pchisq(D, df=df, lower.tail=FALSE))
  D <- list(dispersion=sig2, pvalue=rho, df=df)
  return(D) 
}
