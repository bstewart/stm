topicLasso <- function(formula, data, stmobj=NULL, subset=NULL,
                      omit.var=NULL,
                      family="gaussian", main="Topic Effects on Outcome",
                      xlab=expression("Lower Outcome                  Higher Outcome"),
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
  
  linmod <- cv.glmnet(x=X,y=y, family=family, standardize=standardize, nfolds=nfolds, ...)
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
