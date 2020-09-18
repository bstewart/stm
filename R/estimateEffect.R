#' Estimates regressions using an STM object
#' 
#' Estimates a regression where documents are the units, the outcome is the
#' proportion of each document about a topic in an STM model and the covariates
#' are document-meta data.  This procedure incorporates measurement uncertainty
#' from the STM model using the method of composition.
#' 
#' This function performs a regression where topic-proportions are the outcome
#' variable.  This allows us to conditional expectation of topic prevalence
#' given document characteristics.  Use of the method of composition allows us
#' to incorporate our estimation uncertainty in the dependent variable. Mechanically
#' this means we draw a set of topic proportions from the variational posterior,
#' compute our coefficients, then repeat.  To compute quantities of interest we
#' simulate within each batch of coefficients and then average over all our results.
#' 
#' The formula specifies the nature of the linear model.  On the left hand-side
#' we use a vector of integers to indicate the topics to be included as outcome
#' variables.  If left blank then the default of all topics is used. On the
#' right hand-side we can specify a linear model of covariates including
#' standard transformations.  Thus the model \code{2:4 ~ var1 + s(var2)} would
#' indicate that we want to run three regressions on Topics 2, 3 and 4 with
#' predictor variables \code{var1} and a b-spline transformed \code{var2}.  We
#' encourage the use of spline functions for non-linear transformations of
#' variables.
#' 
#' The function allows the user to specify any variables in the model.
#' However, we caution that for the assumptions of the method of composition to
#' be the most plausible the topic model should contain at least all the
#' covariates contained in the \code{estimateEffect} regression.  However the
#' inverse need not be true.  The function will automatically check whether the
#' covariate matrix is singular which generally results from linearly dependent
#' columns.  Some common causes include a factor variable with an unobserved
#' level, a spline with degrees of freedom that are too high, or a spline with
#' a continuous variable where a gap in the support of the variable results in
#' several empty basis functions.  In these cases the function will still
#' estimate by adding a small ridge penalty to the likelihood.  However, we
#' emphasize that while this will produce an estimate it is only identified by
#' the penalty.  In many cases this will be an indication that the user should
#' specify a different model.
#' 
#' The function can handle factors and numeric variables.  Dates should be
#' converted to numeric variables before analysis.
#' 
#' We offer several different methods of incorporating uncertainty.  Ideally we
#' would want to use the covariance matrix that governs the variational
#' posterior for each document (\eqn{\nu}).  The updates for the global
#' parameters rely only on the sum of these matrices and so we do not store
#' copies for each individual document.  The default uncertainty method
#' \code{Global} uses an approximation to the average covariance matrix formed
#' using the global parameters.  The uncertainty method \code{Local} steps
#' through each document and updates the parameters calculating and then saving
#' the local covariance matrix.  The option \code{None} simply uses the map
#' estimates for \eqn{\theta} and does not incorporate any uncertainty.  We
#' strongly recommend the \code{Global} approximation as it provides the best
#' tradeoff of accuracy and computational tractability.
#' 
#' Effects are plotted based on the results of \code{\link{estimateEffect}}
#' which contains information on how the estimates are constructed.  Note that
#' in some circumstances the expected value of a topic proportion given a
#' covariate level can be above 1 or below 0.  This is because we use a Normal
#' distribution rather than something constrained to the range between 0 and 1.
#' If a continuous variable goes above 0 or 1 within the range of the data it
#' may indicate that a more flexible non-linear specification is needed (such
#' as using a spline or a spline with greater degrees of freedom).
#' 
#' @param formula A formula for the regression.  It should have an integer or
#' vector of numbers on the left-hand side and an equation with covariates on
#' the right hand side.  See Details for more information.
#' @param stmobj Model output from STM
#' @param metadata A dataframe where all predictor variables in the formula can
#' be found. If \code{NULL} R will look for the variables in the global
#' namespace.  It will not look for them in the \code{STM} object which for
#' memory efficiency only stores the transformed design matrix and thus will
#' not in general have the original covariates.
#' @param uncertainty Which procedure should be used to approximate the
#' measurement uncertainty in the topic proportions.  See details for more
#' information.  Defaults to the Global approximation.
#' @param documents If uncertainty is set to \code{Local}, the user needs to
#' provide the documents object (see \code{\link{stm}} for format).
#' @param nsims The number of simulated draws from the variational posterior.
#' Defaults to 25.  This can often go even lower without affecting the results
#' too dramatically.
#' @param prior This argument allows the user to specify a ridge penalty to be
#' added to the least squares solution for the purposes of numerical stability.
#' If its a scalar it is added to all coefficients.  If its a matrix it should
#' be the prior precision matrix (a diagonal matrix of the same dimension as
#' the \code{ncol(X)}).  When the design matrix is collinear but this argument
#' is not specified, a warning will pop up and the function will estimate with
#' a small default penalty.
#' @return \item{parameters}{A list of K elements each corresponding to a
#' topic.  Each element is itself a list of n elements one per simulation.
#' Each simulation contains the MLE of the parameter vector and the variance
#' covariance matrix} \item{topics}{The topic vector} \item{call}{The original
#' call} \item{uncertainty}{The user choice of uncertainty measure}
#' \item{formula}{The formula object} \item{data}{The original user provided
#' meta data.} \item{modelframe}{The model frame created from the formula and
#' data} \item{varlist}{A variable list useful for mapping terms with columns
#' in the design matrix}
#' @seealso \code{\link{plot.estimateEffect}} \code{\link{summary.estimateEffect}}
#' @examples
#' 
#' #Just one topic (note we need c() to indicate it is a vector)
#' prep <- estimateEffect(c(1) ~ treatment, gadarianFit, gadarian)
#' summary(prep)
#' plot(prep, "treatment", model=gadarianFit, method="pointestimate")
#' 
#' #three topics at once
#' prep <- estimateEffect(1:3 ~ treatment, gadarianFit, gadarian)
#' summary(prep)
#' plot(prep, "treatment", model=gadarianFit, method="pointestimate")
#' 
#' #with interactions
#' prep <- estimateEffect(1 ~ treatment*s(pid_rep), gadarianFit, gadarian)
#' summary(prep)
#' @export
estimateEffect <- function(formula,
                     stmobj, metadata=NULL,
                     uncertainty=c("Global", "Local", "None"), documents=NULL,
                     nsims=25, prior=NULL) {
  origcall <- match.call()
  thetatype <- match.arg(uncertainty)
  if(thetatype=="None") nsims <- 1 #override nsims for no uncertainty
  
  if(!is.null(documents)) {
    # Convert the corpus to the internal STM format
    args <- asSTMCorpus(documents, data=metadata)
    documents <- args$documents
    metadata <- args$data
  }
  
  ##
  #Step 1: Extract the formula and do some error checking
  ##
  if(!inherits(formula,"formula")) stop("formula must be a formula object.")
  if(!is.null(metadata) & !is.data.frame(metadata)) metadata <- as.data.frame(metadata)
  termobj <- terms(formula, data=metadata)
  if(attr(termobj, "response")==1){
    #if a response is specified we have to parse it and remove it.
    # as.character of a formula turns
    # dv ~ iv
    # into:  c("~", "dv", "iv")
    response <- as.character(formula)[2] #second object is the response in this cases
    K <- eval(parse(text=response))
    if(!(posint(K) && max(K)<=stmobj$settings$dim$K)) stop("Topics specified as response in formula must be a set of positive integers equal to or less than the number of topics in the model.")   
    #now we reconstruct the formula removing the response
    formula <- formula(paste(as.character(formula)[c(1,3)], collapse = " "))
      #the above used to be the below code but the use got deprecated.
      #as.formula(as.character(formula)[c(1,3)])
    termobj <- terms(formula, data=metadata)
  } else {
    K <- 1:stmobj$settings$dim$K
  }
  mf <- model.frame(termobj, data=metadata)
  xmat <- model.matrix(termobj,data=metadata)
  varlist <- all.vars(termobj)
  if(!is.null(metadata)) {
    data <- metadata[, varlist, drop=FALSE]
  } else {
    templist <- list()
    for(i in 1:length(varlist)) {
      templist[[i]] <- get(varlist[i])
    }
    data <- data.frame(templist)
    names(data) <- varlist
    rm(templist)
  }
  metadata <- data
  rm(data)
  ##
  #Step 2: Compute the QR decomposition
  ##
  # all the models here are essentially just OLS regressions
  # becuase we will run them many times we want to cache the 
  # expensive components in advance.
  if(!is.null(prior)) {
    if(!is.matrix(prior)) {
      prior <- diag(prior, nrow=ncol(xmat))
    } 
    if(ncol(prior)!=ncol(xmat)) stop("number of columns in prior does not match columns in design matrix")
    prior.pseudo <- chol(prior)
    xmat <- rbind(xmat,prior.pseudo)
  }
  qx <- qr(xmat)
  if(qx$rank < ncol(xmat)) {
    prior <- diag(1e-5, nrow=ncol(xmat))
    prior.pseudo <- chol(prior)
    xmat <- rbind(xmat,prior.pseudo)
    qx <- qr(xmat)
    warning("Covariate matrix is singular.  See the details of ?estimateEffect() for some common causes.
             Adding a small prior 1e-5 for numerical stability.")
  }
  ##  
  #Step 3: Calculate Coefficients
  ##

  storage <- vector(mode="list", length=length(K))
  for(i in 1:nsims) {
    # 3a) simulate theta
    if(thetatype=="None") thetasims <- stmobj$theta
    else {
      thetasims <- thetaPosterior(stmobj, nsims=1, type=thetatype, documents=documents)
      thetasims <- do.call(rbind, thetasims)
    }
    # 3b) calculate model
    for(k in K) {
      #lm.mod <- lm(thetasims[,k]~ xmat -1)
      #storage[[which(k==K)]][[i]] <- list(coef=coef(lm.mod),vcov=vcov(lm.mod))
      lm.mod <- qr.lm(thetasims[,k], qx)
      storage[[which(k==K)]][[i]] <- summary.qr.lm(lm.mod)      
    }
  }
    
  ##
  #Step 4: Return Values
  ##
  # 4a) Give it an S3 class
  # 4b) Return the terms object, the call, the formula, the data etc.
  # 4c) anything else we need for a summary() type function
  toreturn <- list(parameters=storage, topics=K,
                   call=origcall, uncertainty=thetatype, formula=formula, data=metadata,
                   modelframe=mf, varlist=varlist)
  class(toreturn) <- "estimateEffect"
  return(toreturn)
}

# A function for performing simple linear regression with a cached QR decomposition
# this should be lighter weight than lm().  Works with summary.qr.lm() to give
# vcov calculations etc.
qr.lm <- function(y, qx) {
  if(length(y)!=nrow(qx$qr)) {
    #probably don't match because of a prior
    if(length(y)!=(nrow(qx$qr)-ncol(qx$qr))) stop("number of covariate observations does not match number of docs")
    #if it passes this check its the prior. thus
    y <- c(y,rep(0, ncol(qx$qr)))
  }
  beta <- solve.qr(qx, y)
  residuals <- qr.resid(qx,y)
  fitted.values <- qr.fitted(qx,y)
  df.residual <- length(fitted.values) - qx$rank
  out <- list(coefficients=beta, residuals=residuals, 
              fitted.values=fitted.values, 
              df.residual=df.residual, rank=qx$rank, qr=qx)
  out 
}
#this function rewrites the summary.lm() function
# to calculate from our reduced regression
summary.qr.lm <- function (object) {
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  
  Qr <- object$qr
  n <- nrow(Qr$qr)
  p1 <- 1L:p
  r <- z$residuals
  f <- z$fitted.values
  
  mss <- ifelse(attr(z$terms, "intercept"), sum((f - mean(f))^2), sum(f^2)) 
  rss <- sum(r^2)
  
  resvar <- rss/rdf
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar)
  est <- z$coefficients[Qr$pivot[p1]]
  sigma <- sqrt(resvar)
  list(est=est, vcov=(sigma^2 * R))
}


#'Summary for estimateEffect
#'
#'Create a summary regression table similar to those produced for \code{lm}
#'
#'This function along with \code{\link{print.summary.estimateEffect}} creates
#'regression tables that look like typically summaries you see in R.  In general
#'we recommend that you use non-linearities such as splines via function like
#'\code{\link{s}} and in those circumstances the tables are not particularly
#'interpretable.  
#'
#'Confidence intervals are calculated by using draws from the covariance matrix
#'of each simulation to estimate the standard error.  Then a t-distribution approximation
#'is applied to calculate the various quantities of interest.
#'
#'@param object an object of class \code{"estimateEffect"}, usually a result of a call to
#'\code{\link{estimateEffect}}
#'@param topics a vector containing the topic numbers for each a summary is to be calculated.
#'Must be contained in the original \code{estimateEffect} object
#'@param nsim the number of simulations to use per parameter set to calculate the standard error.
#'Defaults to 500
#'@param ... further arguments passed to or from other methods
#'
#'@seealso \code{\link{estimateEffect}} \code{\link{plot.estimateEffect}}
#'@method summary estimateEffect  
#'@aliases summary.estimateEffect print.summary.estimateEffect
#'@export
summary.estimateEffect <- function(object, topics=NULL, nsim=500, ...) {
  if(is.null(topics)) topics <- object$topics
  if(any(!(topics %in% object$topics))) {
    stop("Some topics specified with the topics argument are not available in this estimateEffect object.")
  }
  tables <- vector(mode="list", length=length(topics))
  for(i in 1:length(topics)) {
    topic <- topics[i]
    sims <- lapply(object$parameters[[which(object$topics==topic)]], function(x) rmvnorm(nsim, x$est, x$vcov))
    sims <- do.call(rbind,sims)
    est<- colMeans(sims)
    se <- sqrt(apply(sims,2, stats::var))
    tval <- est/se
    rdf <- nrow(object$data) - length(est)
    p <- 2 * stats::pt(abs(tval), rdf, lower.tail = FALSE)
    
    coefficients <- cbind(est, se, tval, p)
    rownames(coefficients) <- attr(object$parameters[[1]][[1]]$est, "names") 
    colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    tables[[i]] <- coefficients
  }
  out <- list(call=object$call, topics=topics, tables=tables)
  class(out) <- "summary.estimateEffect"
  return(out)
}

#'@method print summary.estimateEffect  
#'@export
print.summary.estimateEffect <- function(x, digits = max(3L, getOption("digits") - 3L), 
                                         signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  for(i in 1:length(x$tables)) {
    cat(sprintf("\nTopic %i:\n", x$topics[i]))
    cat("\nCoefficients:\n")
    coefs <- x$tables[[i]]
    stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    cat("\n")
  }
  invisible(x)
}