## Structural Topic Model
# this is a wrapper around workhorse function stm.control() 
# the primary purpose here is to format arguments and do some light preprocessing.

stm <- function(documents, vocab, K, 
                prevalence, content, data=NULL,
                init.type=c("LDA", "DMR","Random"), seed=NULL, 
                max.em.its=100, emtol=1e-5,
                verbose=TRUE, reportevery=5, keepHistory=FALSE,  
                LDAbeta=TRUE, interactions=TRUE,
                gamma.prior=c("Pooled", "L1"), sigma.prior=0,
                kappa.prior=c("Jeffreys", "L1"), control=list())  {
  
  #Match Arguments and save the call
  init.type <- match.arg(init.type)
  Call <- match.call()
  
  #Documents
  if(missing(documents)) stop("Must include documents")
  if(!is.list(documents)) stop("documents must be a list, see documentation.")
  if(!all(unlist(lapply(documents, is.matrix)))) stop("Each list element in documents must be a matrix. See documentation.")
  
  N <- length(documents)
  
  #Extract and Check the Word indices
  wcounts <- doc.to.ijv(documents)
  wcounts <- aggregate(wcounts$v, by=list(wcounts$j), FUN=sum)
  V <- length(wcounts$Group.1)  
  if(!posint(wcounts$Group.1)) {
    stop("Word indices are not positive integers")
  } 
  if(!isTRUE(all.equal(wcounts$Group.1,1:V))) {
    stop("Word indices must be sequential integers starting with 1.")
  } 
  
  #Check the Vocab vector against the observed word indices
  if(length(vocab)!=V) stop("Vocab length does not match observed word indices")
  
  #Check the Number of Topics
  if(missing(K)) stop("K, the number of topics, is required.")
  if(!(posint(K) && length(K)==1 && K>1)) stop("K must be a positive integer greater than 1.")
  if(K==2) warning("K=2 is equivalent to a unidimensional scaling model which you may prefer.")
  
  #Iterations, Verbose etc.
  if(!(length(max.em.its)==1 & posint(max.em.its))) stop("Max EM iterations must be a single positive integer")
  if(!is.logical(verbose)) stop("verbose must be a logical.")
  
  ##
  # A Function for processing prevalence-covariate design matrices
  makeTopMatrix <- function(x, data=NULL) {
    #is it a formula?
    if(inherits(x,"formula")) {
      termobj <- terms(x, data=data)
      if(attr(termobj, "response")==1) stop("Response variables should not be included in prevalence formula.")
      xmat <- model.matrix(termobj,data=data)
      return(xmat)
    }
    if(is.matrix(x)) {
      #Does it have an intercept in first column?
      if(isTRUE(all.equal(x[,1],rep(1,nrow(x))))) return(x) 
      else return(cbind(1,x))
    }
  }
  
  ###
  #Now we parse both sets of covariates
  ###
  if(!missing(prevalence)) {
    if(!is.matrix(prevalence) & !inherits(prevalence, "formula")) stop("Prevalence Covariates must be specified as a model matrix or as a formula")
    xmat <- makeTopMatrix(prevalence,data)
    if(nrow(na.omit(xmat)) != length(documents)) stop("Complete cases in prevalence covariate does not match the number of documents.")
  } else {
    xmat <- NULL
  }
  
  if(!missing(content)) {
    if(inherits(content, "formula")) {
      termobj <- terms(content, data=data)
      if(attr(termobj, "response")==1) stop("Response variables should not be included in content formula.")
      if(nrow(attr(termobj, "factors"))!=1) stop("Currently content can only contain one variable.")
      if(is.null(data)) {
        yvar <- eval(attr(termobj, "variables"))[[1]]
      } else {
        char <- rownames(attr(termobj, "factors"))[1]
        yvar <- data[[char]]
      }
      yvar <- as.factor(yvar)
    } else {
      yvar <- as.factor(content)
    }
      yvarlevels <- levels(yvar)
      betaindex <- as.numeric(yvar)
  } else{
    yvarlevels <- NULL
    betaindex <- rep(1, length(documents))
  }
  A <- length(unique(betaindex)) #define the number of aspects
  
  #Checks for Dimension agreement
  ny <- length(betaindex)
  nx <- ifelse(is.null(xmat), N, nrow(xmat))
  if(N!=nx | N!=ny) stop(paste("number of observations in content covariate (",ny,
                                  ") prevalence covariate (",
                                  nx,") and documents (",N,") are not all equal.",sep=""))
  
  #Some additional sanity checks
  if(!is.logical(LDAbeta)) stop("LDAbeta must be logical")
  if(!is.logical(interactions)) stop("Interactions variable must be logical")
  if(sigma.prior < 0 | sigma.prior > 1) stop("sigma.prior must be between 0 and 1")

  ###
  # Now Construct the Settings File
  ###
  settings <- list(dim=list(K=K, A=A, 
                            V=V, N=N, wcounts=wcounts),
                   verbose=verbose,
                   topicreportevery=reportevery,
                   keepHistory=keepHistory,
                   convergence=list(max.em.its=max.em.its, em.converge.thresh=emtol, topwords=20, topwordits=0),
                   covariates=list(X=xmat, betaindex=betaindex, yvarlevels=yvarlevels),
                   gamma=list(mode=match.arg(gamma.prior), prior=NULL, enet=1),
                   sigma=list(prior=sigma.prior),
                   kappa=list(LDAbeta=LDAbeta, interactions=interactions, 
                              fixedintercept=TRUE, mstep=list(maxit=3, tol=.99)),
                   tau=list(mode=match.arg(kappa.prior), prior=NULL, tol=1e-5,
                            enet=1,nlambda=500, lambda.min.ratio=.0001, ic.k=2),
                   init=list(mode=init.type, 
                             userinit=NULL), 
                   seed=seed)
  if(settings$gamma$mode=="L1") {
    if(!require(glmnet) | !require(Matrix)) stop("To use L1 penalization please install glmnet and Matrix")
    if(ncol(xmat)<=2) stop("Cannot use L1 penalization in prevalence model with 2 or fewer covariates.")
  }
  if(settings$tau$mode=="L1") {
    if(!require(glmnet) | !require(Matrix)) stop("To use L1 penalization please install glmnet and Matrix")    
    settings$tau$maxit <- 1e8
    settings$tau$tol <- 1e-6
  }
  
  ###
  # Fill in some implied arguments.
  ###
  
  #Is there a covariate on top?
  if(missing(prevalence)) {
    settings$gamma$mode <- "CTM" #without covariates has to be estimating the mean.
  } 
  
  #Is there a covariate on the bottom?
  if(missing(content)) {
    settings$kappa$interactions <- FALSE #can't have interactions without a covariate.
  } else {
    settings$kappa$LDAbeta <- FALSE #can't do LDA topics with a covariate 
  }
  
  ###
  # process arguments in control
  ###
  
  #Full List of legal extra arguments
  legalargs <-  c("tau.maxit", "tau.tol","kappa.mstepmaxit", "kappa.msteptol", 
                  "wordconverge.num", "wordconverge.its", "fixedintercept",
                  "kappa.enet", "nlambda", "lambda.min.ratio", "ic.k", "gamma.enet")
  if (length(control)) {
    indx <- pmatch(names(control), legalargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument %s not matched", names(control)[indx==0L]),
           domain = NA)
    fullnames <- legalargs[indx]
    for(i in fullnames) {
      if(i=="tau.maxit") settings$tau$maxit <- control[[i]]
      if(i=="tau.tol") settings$tau$tol <- control[[i]]
      if(i=="kappa.mstepmaxit") settings$kappa$mstep$maxit <- control[[i]] 
      if(i=="kappa.msteptol") settings$kappa$mstep$tol <- control[[i]] 
      if(i=="wordconverge.num") settings$convergence$topwords <- control[[i]]
      if(i=="wordconverge.its") settings$convergence$topwordits <- control[[i]]
      if(i=="fixedintercept")settings$kappa$fixedintercept <- control[[i]]
      if(i=="kappa.enet") settings$tau$enet <- control[[i]]
      if(i=="nlambda") settings$tau$nlambda <- control[[i]]
      if(i=="lambda.min.ratio") settings$tau$lambda.min.ratio <- control[[i]]
      if(i=="ic.k") settings$tau$ic.k <- control[[i]]
      if(i=="gamma.enet") settings$gamma$enet <- control[[i]]
    }
  }
  
  ###
  # Process the Seed
  ###
  if(is.null(settings$seed)) {
    #if there is no seed, choose one and set it, recording for later
    seed <- floor(runif(1)*1e7) 
    set.seed(seed)
    settings$seed <- seed
  } else {
    #otherwise just use the provided seed.
    set.seed(settings$seed)
  }
  
  settings$call <- Call
  ###
  # Finally run the actual model
  ###
  return(stm.control(documents, vocab, settings))
}
