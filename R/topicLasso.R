#' Plot predictions using topics
#' 
#' Use the \pkg{glmnet} package to plot LASSO based estimates of relationship
#' between an arbitrary dependent variable with topics and additional variables
#' as predictors.  This function is experimental (see below).
#' 
#' This function is used for estimating the most important topical predictors
#' of an arbitrary outcome.  The idea is to run an L1 regularized regression
#' using \code{\link{cv.glmnet}} in the \pkg{glmnet} package where the
#' document-level dependent variable is chosen by the user and the predictors
#' are the document-topic proportions in the \code{\link{stm}} model along with
#' any other variables of interest.
#' 
#' The function uses cross-validation to choose the regularization parameter
#' and generates a plot of which loadings were the most influential in
#' predicting the outcome.  It also invisibly returns the glmnet model so that
#' it can be used for prediction.
#' 
#' NOTE: This function is still very experimental and may have stability
#' issues.  If stability issues are encountered see the documentation in
#' \pkg{glmnet} for arguments that can be passed to improve convergence.  Also,
#' it is unlikely to work well with multivariate gaussian or multinomial
#' families.
#' 
#' @param formula Formula specifying the dependent variable and additional
#' variables to included in the LASSO beyond the topics present in the stmobj.
#' Just pass a 1 on the right-hand side in order to run without additional
#' controls.
#' @param data Data file containing the dependent variable. Typically will be
#' the metadata file used in the stm analysis. It must have a number of rows
#' equal to the number of documents in the stmobj.
#' @param stmobj The STM object, and output from the \code{stm} function.
#' @param subset A logical statement that will be used to subset the corpus.
#' @param omit.var Pass a character vector of variable names to be excluded
#' from the plot.  Note this does not exclude them from the calculation, only
#' the plot.
#' @param family The family parameter used in \code{\link{glmnet}}.  See
#' explanation there.  Defaults to "gaussian"
#' @param main Character string for the main title.
#' @param xlab Character string giving an x-axis label.
#' @param labeltype Type of example words to use in labeling each topic. See
#' \code{\link{labelTopics}}. Defaults to "prob".
#' @param seed The random seed for replication of the cross-validation samples.
#' @param xlim Width of the x-axis.
#' @param standardize Whether to standardize variables. Default is FALSE, which
#' is different from the \pkg{glmnet} default because the topics are already
#' standardized.  Note that glmnet standardizes the variables by default but
#' then projects them back to their original scales before reporting
#' coefficients.
#' @param nfolds the number of cross-validation folds.  Defaults to 20.
#' @param ...  Additional arguments to be passed to glmnet.  This can be useful
#' for addressing convergence problems.
#' @seealso \pkg{glmnet}
#' @references Friedman, Jerome, Trevor Hastie, and Rob Tibshirani.
#' "Regularization paths for generalized linear models via coordinate descent."
#' Journal of statistical software 33.1 (2010): 1.
#' @examples
#' 
#' \donttest{
#' 
#' #Load the poliblog data
#' data(poliblog5k)
#' #estimate a model with 50 topics
#' stm1 <- stm(poliblog5k.docs, poliblog5k.voc, 50,
#'             prevalence=~rating + blog, data=poliblog5k.meta,
#'             init.type="Spectral")
#' #make a plot of the topics most predictive of "rating"
#' out <- topicLasso(rating ~ 1, family="binomial", data=poliblog5k.meta,stmobj=stm1)
#' #generate some in-sample predictions
#' pred <- predict(out, newx=stm1$theta,type="class")
#' #check the accuracy of the predictions
#' table(pred, poliblog5k.meta$rating)
#' }
#' @export
topicLasso <- function(formula, data, stmobj=NULL, subset=NULL,
                      omit.var=NULL, family="gaussian", 
                      main="Topic Effects on Outcome", 
                      xlab=expression("Lower Outcome Higher Outcome"),
                      labeltype=c("prob", "frex", "lift", "score"),seed=02138,
                      xlim=c(-4,4), standardize=FALSE, nfolds=20, ...) {

  X <- model.matrix(formula, data)
  pvarnames <- colnames(X)
  p <- ncol(X)
  z <- stmobj
  X <- cbind(X,z$theta) 
  #if there is an intercept remove it.
  if(colnames(X)[1]=="(Intercept)") {
    p <- p - 1
    X <- X[,-1]
    pvarnames <- pvarnames[-1]
  }
  mframe <- model.frame(formula,data)
  y <- model.response(mframe)
  
  if(!is.null(subset)) {
    X <- X[subset,]
    y <- y[subset]
  }
  X <- Matrix(X)
  
  
  
  if(!is.null(seed)) set.seed(seed)

  
  #pick the label type- if its a conten covariate fix it
  labeltype <- match.arg(labeltype)
  if(!z$settings$kappa$LDAbeta) labeltype <- "topics"
  
  topiclabs <- labelTopics(z)
  topiclabs <- apply(topiclabs[[labeltype]], 1, commas)
  topiclabs <- sprintf("Topic %i: %s", 1:length(topiclabs), topiclabs)
  
  varnames <- c(pvarnames, topiclabs)
  
  linmod <- glmnet::cv.glmnet(x=X,y=y, family=family, standardize=standardize, nfolds=nfolds, ...)
  loadings <- coef(linmod)[-1]
  
  posind <- which(loadings>0)
  if(length(posind!=0)) {
      cat(c("Positive Effect \n",sprintf("%s \n", 
                                        varnames[posind])))
  } else {
    cat("No Positive Effects. \n")
  }
  negind <- which(loadings<0)
  
  if(length(negind!=0)) {
    cat(c("Negative Effect \n",sprintf("%s \n", 
                                       varnames[negind])))
  } else {
    cat("No Negative Effects. \n")
  }
    
  coefind <- which(loadings!=0)
  drop <- c()
  if(!is.null(omit.var)) {
    for(i in 1:length(omit.var)) {
      drop <- c(drop,agrep(omit.var[i], varnames))
    }
    coefind <- coefind[!(coefind %in% drop)]
  }
  if(length(coefind)!=0) {
    coefval <- loadings[coefind]
    labels <- varnames[coefind]
    
    labels <- labels[order(coefval)]
    values <- sort(coefval)
    plot(x=values, y=1:length(values), pch=19, main=main,
         xlim=xlim, ylab="", type="p", xlab=xlab,
         ylim=c(-1, length(values)+1),yaxt="n", bty="n")
    text(values, 1:length(values), labels, pos=2)
    segments(0, -1, 0, length(values)+1,lty=2)
  } else {
    cat("No Non Zero Loadings")
  }
  return(invisible(linmod))
}
