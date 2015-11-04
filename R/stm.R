## Structural Topic Model
# wrapper for the stm.control function.  

stm <- function(documents, vocab, K, 
                prevalence, content, data=NULL,
                init.type=c("LDA", "Random", "Spectral"), seed=NULL, 
                max.em.its=500, emtol=1e-5,
                verbose=TRUE, reportevery=5,   
                LDAbeta=TRUE, interactions=TRUE, 
                ngroups=1, model=NULL,
                gamma.prior=c("Pooled", "L1"), sigma.prior=0,
                kappa.prior=c("L1", "Jeffreys"), control=list())  {
  
  #Match Arguments and save the call
  init.type <- match.arg(init.type)
  Call <- match.call()
  
  #Documents
  if(missing(documents)) stop("Must include documents")
  if(!is.list(documents)) stop("documents must be a list, see documentation.")
  if(!all(unlist(lapply(documents, is.matrix)))) stop("Each list element in documents must be a matrix. See documentation.")
  if(any(unlist(lapply(documents, function(x) anyDuplicated(x[1,]))))) {
    stop("Duplicate term indices within a document.  See documentation for proper format.")
  }
  N <- length(documents)
  
  #Extract and Check the Word indices
  wcountvec <- unlist(lapply(documents, function(x) rep(x[1,], times=x[2,])),use.names=FALSE)
  #to make this backward compatible we reformulate to old structure.
  wcounts <- list(Group.1=sort(unique(wcountvec)))
  V <- length(wcounts$Group.1)  
  if(!posint(wcounts$Group.1)) {
    stop("Word indices are not positive integers")
  } 
  if(!isTRUE(all.equal(wcounts$Group.1,1:V))) {
    stop("Word indices must be sequential integers starting with 1.")
  } 
  #note we only do the tabulation after making sure it will actually work.
  wcounts$x <- tabulate(wcountvec)
  rm(wcountvec)
  
  #Check the Vocab vector against the observed word indices
  if(length(vocab)!=V) stop("Vocab length does not match observed word indices")
  
  #Check the Number of Topics
  if(missing(K)) stop("K, the number of topics, is required.")
  if(K!=0) {
    #this is the old set of checks
    if(!(posint(K) && length(K)==1 && K>1)) stop("K must be a positive integer greater than 1.")
    if(K==2) warning("K=2 is equivalent to a unidimensional scaling model which you may prefer.")
  } else {
    #this is the special set of checks for Lee and Mimno
    if(init.type!="Spectral") stop("Topic selection method can only be used with init.type='Spectral'")
  }
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
      xmat <- try(sparse.model.matrix(termobj,data=data),silent=TRUE)
      if(class(xmat)=="try-error") stop("Error creating model matrix.
                                        This could be caused by many things including
                                        explicit calls to a namespace within the formula.
                                        Try a simpler formula.")
      propSparse <- 1 - nnzero(xmat)/length(xmat) 
      #if its less than 50% sparse or there are fewer than 50 columns, just convert to a standard matrix
      if(propSparse < .5 | ncol(xmat) < 50) {
        xmat <- as.matrix(xmat)
      }
      return(xmat)
    }
    if(is.matrix(x)) {
      #Does it have an intercept in first column?
      if(isTRUE(all.equal(x[,1],rep(1,nrow(x))))) return(Matrix(x)) 
      else return(cbind(1,Matrix(x)))
    }
  }
  
  ###
  #Now we parse both sets of covariates
  ###
  if(!missing(prevalence)) {
    if(!is.matrix(prevalence) & !inherits(prevalence, "formula")) stop("Prevalence Covariates must be specified as a model matrix or as a formula")
    xmat <- makeTopMatrix(prevalence,data)
    if(is.na(nnzero(xmat))) stop("Missing values in prevalence covariates.")
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
  if(!is.null(model)) {
    if(max.em.its <= model$convergence$its) stop("when restarting a model, max.em.its represents the total iterations of the model 
                                                 and thus must be greater than the length of the original run")
  }
  ###
  # Now Construct the Settings File
  ###
  settings <- list(dim=list(K=K, A=A, 
                            V=V, N=N, wcounts=wcounts),
                   verbose=verbose,
                   topicreportevery=reportevery,
                   convergence=list(max.em.its=max.em.its, em.converge.thresh=emtol),
                   covariates=list(X=xmat, betaindex=betaindex, yvarlevels=yvarlevels),
                   gamma=list(mode=match.arg(gamma.prior), prior=NULL, enet=1),
                   sigma=list(prior=sigma.prior),
                   kappa=list(LDAbeta=LDAbeta, interactions=interactions, 
                              fixedintercept=TRUE, mstep=list(tol=.001, maxit=3),
                              contrast=FALSE),
                   tau=list(mode=match.arg(kappa.prior), tol=1e-5,
                            enet=1,nlambda=250, lambda.min.ratio=.001, ic.k=2,
                            maxit=1e4),
                   init=list(mode=init.type, nits=50, burnin=25, alpha=(50/K), eta=.01,
                             s=.05, p=3000, d.group.size=2000), 
                   seed=seed,
                   ngroups=ngroups)
  if(settings$gamma$mode=="L1") {
    #if(!require(glmnet) | !require(Matrix)) stop("To use L1 penalization please install glmnet and Matrix")
    if(ncol(xmat)<=2) stop("Cannot use L1 penalization in prevalence model with 2 or fewer covariates.")
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
  legalargs <-  c("tau.maxit", "tau.tol", 
                  "fixedintercept","kappa.mstepmaxit", "kappa.msteptol", 
                  "kappa.enet", "nlambda", "lambda.min.ratio", "ic.k", "gamma.enet",
                  "nits", "burnin", "alpha", "eta", "contrast",
                  "rp.s", "rp.p", "rp.d.group.size", "SpectralRP")
  if (length(control)) {
    indx <- pmatch(names(control), legalargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument %s not matched", names(control)[indx==0L]),
           domain = NA)
    fullnames <- legalargs[indx]
    for(i in fullnames) {
      if(i=="tau.maxit") settings$tau$maxit <- control[[i]]
      if(i=="tau.tol") settings$tau$tol <- control[[i]]
      if(i=="fixedintercept")settings$kappa$fixedintercept <- control[[i]]
      if(i=="kappa.enet") settings$tau$enet <- control[[i]]
      if(i=="kappa.mstepmaxit") settings$kappa$mstep$maxit <- control[[i]] 
      if(i=="kappa.msteptol") settings$kappa$mstep$tol <- control[[i]] 
      if(i=="nlambda") settings$tau$nlambda <- control[[i]]
      if(i=="lambda.min.ratio") settings$tau$lambda.min.ratio <- control[[i]]
      if(i=="ic.k") settings$tau$ic.k <- control[[i]]
      if(i=="gamma.enet") settings$gamma$enet <- control[[i]]
      if(i=="nits") settings$init$nits <- control[[i]]
      if(i=="burnin") settings$init$burnin <- control[[i]]
      if(i=="alpha") settings$init$alpha <- control[[i]]
      if(i=="eta") settings$init$eta <- control[[i]]
      if(i=="contrast") settings$kappa$contrast <- control[[i]]
      if(i=="rp.s")  settings$init$s <- control[[i]]
      if(i=="rp.p")  settings$init$p <- control[[i]]
      if(i=="rp.d.group.size")  settings$init$d.group.size <- control[[i]]
      if(i=="SpectralRP" & control[[i]]) settings$init$mode <- "SpectralRP" #override to allow spectral rp mode
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
  return(stm.control(documents, vocab, settings,model))
}
