## Structural Topic Model
# this is a wrapper around workhorse function stm.control() 
# the primary purpose here is to format arguments and do some light preprocessing.
stm <- function(documents, vocab, K, 
                prevalence, content, data=NULL,
                init.type=c("LDA", "DMR","Random"), seed=NULL, 
                max.em.its=100, emtol=1e-5,
                verbose=TRUE, ...) {
  
  #Very First Argument Parses
  init.type <- match.arg(init.type)
  Call <- match.call()
  
  #Documents
  if(missing(documents)) stop("Must include documents")
  if(!is.list(documents)) stop("documents must be a list, see documentation.")
  if(!all(unlist(lapply(documents, is.matrix)))) stop("Each list element in documents must be a matrix. See documentation.")
  N <- length(documents)
  Vindex <- sort(unique(unlist(lapply(documents, function(x) x[1,]))))
  V <- length(Vindex)
  if(!posint(Vindex)) stop("Word indices are not positive integers")
  if(!isTRUE(all.equal(Vindex,1:V))) stop("Word indices must be sequential integers starting with 1.")
  
  #Vocab
  if(length(vocab)!=V) stop("Vocab length does not match observed word indices")
  
  #K
  if(missing(K)) stop("K, the number of topics, is required.")
  if(!(posint(K) && length(K)==1 && K>1)) stop("K must be a positive integer greater than 1.")
  if(K==2) warning("K=2 is equivalent to a unidimensional scaling model which you may prefer.")
  
  #Iterations, Verbose etc.
  if(!(length(max.em.its)==1 & posint(max.em.its))) stop("Max EM iterations must be a single positive integer")
  if(!is.logical(verbose)) stop("verbose must be a logical.")
  
  ####
  # Function for processing Top-Covariates.
  ####
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
  #Now we Parse Covariates on the Top and Bottom
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
  
  #Checks for Dimension agreement
  ny <- length(betaindex)
  nx <- ifelse(is.null(xmat), N, nrow(xmat))
  if(N!=nx | N!=ny) stop(paste("number of observations in content covariate (",ny,
                                  ") prevalence covariate (",
                                  nx,") and documents (",N,") are not all equal.",sep=""))
  
  ###
  # Now Construct the Settings File
  ###
  wcounts <- doc.to.ijv(documents)
  wcounts <- aggregate(wcounts$v, by=list(wcounts$j), FUN=sum)
  settings <- list(dim=list(K=K, A=length(unique(betaindex)), V=V, N=N, wcounts=wcounts),
                   verbose=verbose,
                   topicreportevery=5,
                   keepHistory=FALSE,
                   convergence=list(max.em.its=max.em.its, em.converge.thresh=emtol, topwords=20, topwordits=0),
                   covariates=list(X=xmat, betaindex=betaindex, yvarlevels=yvarlevels),
                   gamma=list(mode="Pooled", prior=NULL),
                   sigma=list(prior=0),
                   kappa=list(LDAbeta=FALSE, interactions=TRUE, 
                              mstep=list(maxit=3, tol=.99)),
                   tau=list(mode="Jeffreys", prior=NULL, maxit=20, tol=1e-5),
                   init=list(mode=init.type, 
                             userinit=NULL), 
                   seed=seed)
  ###
  # Fill in some implied arguments.
  ###
  
  #Is there a covariate on top?
  if(missing(prevalence)) {
    settings$gamma$mode <- "CTM"
  } 
  
  #Is there a covariate on the bottom?
  if(missing(content)) {
    settings$dim$A <- 1 #no covariate means we only have one aspect.
    settings$kappa$LDAbeta <- TRUE #do LDA topics unless instructed otherwise
    settings$kappa$interactions <- FALSE #can't have interactions without a covariate.
  } 
  
  ###
  # Process Extra Arguments in the Dots Expansion
  ###
  #Full List of legal extra arguments
  legalargs <-  c("LDAbeta", "interactions","gammamode", "sigmaprior", 
                  "taumode", "tauprior", "taumaxit", "tautol",
                  "kappamstepmaxit", "kappamsteptol", "reportevery", 
                  "keepHistory", "userinit", "wordconverge.num", "wordconverge.its")
  extraArgs <- list(...)
  if (length(extraArgs)) {
    indx <- pmatch(names(extraArgs), legalargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument %s not matched", names(extraArgs)[indx==0L]),
           domain = NA)
    fullnames <- legalargs[indx]
    for(i in fullnames) {
      if(i=="LDAbeta") {
        LDAbeta <- extraArgs[[i]]
        if(!is.logical(LDAbeta)) stop("LDAbeta must be logical")
        settings$kappa$LDAbeta <- LDAbeta
      }
      if(i=="interactions") {
        interactions <- extraArgs[[i]]
        if(!is.logical(interactions)) stop("Interactions variable must be logical")
        settings$kappa$interactions <- interactions        
      }
      if(i=="gammamode") {
        gammamode <- extraArgs[[i]]
        options <- c("Pooled", "GL")
        userchoice <- match(gammamode, options, nomatch=0)
        if(userchoice==0) stop("unrecognized gammamode")
        if(userchoice==2) {
          if(!require(gamlr)) stop("To use the gamma lasso please install gamlr")
        }
        else settings$gamma$mode <- options[userchoice]
      }
      if(i=="sigmaprior") {
        sigmaprior <- extraArgs[[i]]
        if(sigmaprior >= 0 & sigmaprior <= 1) settings$sigma$prior <- sigmaprior
        else stop("sigmaprior must be between 0 and 1")
      }
      if(i=="taumode") {
        taumode <- extraArgs[[i]]
        options <- c("Pooled", "Fixed", "Jeffreys", "GL")
        userchoice <- match(taumode, options, nomatch=0)
        if(userchoice==0) stop("unrecognized taumode")
        if(userchoice==4) stop("Gamma Lasso option coming soon.  
                               textir package recently changed 
                               their functionality and we haven't 
                               yet updated accordingly")
        settings$tau$mode <- options[userchoice]
      }
      if(i=="tauprior") {
        settings$tau$prior <- extraArgs[[i]]
      }
      if(i=="taumaxit") {
        settings$tau$maxit <- extraArgs[[i]]
      }
      if(i=="tautol") {
        settings$tau$tol <- extraArgs[[i]]
      }
      if(i=="kappamstepmaxit") {
        settings$kappa$mstep$maxit <- extraArgs[[i]] 
      }
      if(i=="kappamsteptol") {
        settings$kappa$mstep$tol <- extraArgs[[i]] 
      }
      if(i=="reportevery") {
        settings$topicreportevery  <- extraArgs[[i]]
      }
      if(i=="keepHistory") {
        settings$keepHistory <- extraArgs[[i]]
      }
      if(i=="userinit") {
        if(init.type=="User") {
          if(inherits(extraArgs[[i]], "STMR")) {
            settings$init$userinit <- extraArgs[[i]]$beta$logbeta
          } else {
            settings$init$userinit <- extraArgs[[i]]
          }
        } else{
          warning("User Initialization provided but initialization type not set to User.  Using method specified by init.type",immediate=TRUE)
        }
      }
      if(i=="wordconverge.num") {
        settings$convergence$topwords <- wordconverge.num
      }
      if(i=="wordconverge.its") {
        settings$convergence$topwordits <- wordconverge.its
      }
    }
  }
  
  ###
  # Deal with the Seed
  ###
  if(is.null(settings$seed)) {
    seed <- floor(runif(1)*1e7) 
    set.seed(seed)
    settings$seed <- seed
  } else {
    set.seed(settings$seed)
  }
  
  settings$call <- Call
  ###
  # Finally run the actual model
  ###
  return(stm.control(documents, vocab, settings))
}
