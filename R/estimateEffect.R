estimateEffect <- function(formula,
                     stmobj, metadata=NULL,
                     uncertainty=c("Global", "Local", "None"), documents=NULL,
                     nsims=25, prior=NULL) {
  origcall <- match.call()
  thetatype <- match.arg(uncertainty)
  if(thetatype=="None") nsims <- 1 #override nsims for no uncertaintys
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
    formula <- as.formula(as.character(formula)[c(1,3)])
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
  
