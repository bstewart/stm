#' Fit New Documents
#' 
#' A function for predicting thetas for an unseen document based on the previously fit model.
#' 
#' Due to the existence of the metadata in the model, this isn't as simple as in models
#' without side information such as Latent Dirichlet Allocation. There are four scenarios:
#' models without covariate information, models with prevalence covariates only, models with
#' content covariates only and models with both.  When there is not covariate information the
#' choice is essentially whether or not to use prior information.
#' 
#' We offer three types of choices (and may offer more in the future):
#' \describe{
#' \item{"None"}{No prior is used.  In the prevalence case this means that the model simply
#' maximizes the likelihood of seeing the words given the word-topic distribution.  This will
#' in general produce more sharply peaked topic distributions than the prior. This can be used
#' even without the covariates.  This is not an option for topical content covariate models.  If
#' you do not observe the topical content covariate, use the "Average" option.} 
#' \item{"Average"}{We use a prior that is based on the average over the documents in the training
#' set.  This does not require the unseen documents to observe the covariates.  In a model that originally
#' had covariates we need to adjust our estimate of the variance-covariance matrix sigma to accomodate that
#' we no longer have the covariate information.  So we recalculate the variance based on what it would have
#' been if we didn't have any covariates.  This helps avoid an edge case where the covariates are extremely
#' influential and we don't want that strength applied to the new covariate-less setting.  In the case of 
#' content covariates this essentially use the \code{\link{sageLabels}} approach to create a
#' marginalized distribution over words for each topic.
#' }
#' \item{"Covariate"}{We use the same covariate driven prior that existed in the original model.
#' This requires that the test covariates be observed for all previously unseen documents.}
#' }
#' 
#' If you fit a document that was used during training with the options to replicate the initial
#' \code{\link{stm}} model fit you will not necessarily get exactly the same result.  \code{\link{stm}}
#' updates the topic-word distributions last so they may shifted since the document-topic proportions
#' were updated.  If the original model converged, they should be very close.
#' 
#' By default the function returns only the MAP estimate of the normalized document-topic proportions
#' theta.  By selecting \code{returnPrior=TRUE} you can get the various model parameters used to complete
#' the fit.  By selecting \code{returnPosterior=TRUE} you can get the full variational posterior.  Please
#' note that the full variational posterior is very memory intensive.  For a sense of scale it requires an
#' extra \eqn{K^2 + K \times (V'+1) + 1} doubles per document where V' is the number of unique tokens in the document. 
#' 
#' \strong{Testing:} Getting the prevalence covariates right in the unseen documents can be tricky.  However
#' as long as you leave \code{test} set to \code{TRUE} the code will automatically run a test to make sure
#' that everything lines up.  See the internal function \code{\link{makeDesignMatrix}} for more on what is 
#' going on here.
#' 
#' \strong{Passing a Design Matrix}  Advanced users may wish to circumvent this process and pass their
#' own design matrix possibly because they used their own function for transforming the original input
#' variables.  This can be done by passing the design matrix using the \code{designMatrix} argument
#' The columns need to match the ordering of the design matrix for the original \code{stm} object.  
#' The design matrix in an stm model called \code{stmobj} can be found in \code{stmobj$settings$covariates$X} 
#' which can in turn be used to check that you have formatted your result correctly. If you are going to 
#' try this we recommend that you read the documentation for \code{\link{makeDesignMatrix}} to understand
#' some of the challenges involved.  
#' 
#' If you want even more fine-grained control we recommend you directly use the 
#' optimization function \code{\link{optimizeDocument}}
#' 
#' @param model the originally fit STM object.
#' @param documents the new documents to be fit. These documents must be in the stm format and
#'  be numbered in the same way as the documents in the original model with the same dimension of vocabulary.
#'  See \code{\link{alignCorpus}} or the \pkg{quanteda} feature \link[quanteda]{dfm_select} 
#'  for a way to do this.  
#' @param newData the metadata for the prevalence prior which goes with the unseen documents. As in
#' the original data this cannot have any missing data.
#' @param origData the original metadata used to fit the STM object. 
#' @param prevalence the original formula passed to prevalence when \code{stm} was called. The function
#' will try to reconstruct this.
#' @param betaIndex a vector which indicates which level of the content covariate is used
#' for each unseen document. If originally passed as a factor, this can be passed as a factor or 
#' character vector as well but it must not have any levels not included in the original factor.
#' @param prevalencePrior three options described in detail below.  Defaults to "Average" when
#' \code{data} is not provided and "Covariate" when it is.
#' @param contentPrior two options described in detail below. Defaults to "Average" when 
#' \code{betaIndex} is not provided and "Covariate" when it is.
#' @param returnPosterior the function always returns the posterior mode of theta
#' (document-topic proportions),  If set to TRUE this will return the full variational
#' posterior.  Note that this will return a dense K-by-K matrix for every document
#' which can be very memory intensive if you are processing a lot of documents.
#' @param returnPriors the function always returns the options that were set for the prior
#'  (either by the user or chosen internally by the defaults).  In the case of content covariates
#'  using the covariate prior this will be a set of indices to the original beta matrix so as
#'  not to make the object too large.
#' @param designMatrix an option for advanced users to pass an already constructed design matrix for
#' prevalence covariates.  This will override the options in \code{newData}, \code{origData} and
#' \code{test}.  See details below- please do not attempt to use without reading carefully.
#' @param test a test of the functions ability to reconstruct the original functions.
#' @param verbose Should a dot be printed every time 1 percent of the documents are fit.
#' 
#' @return an object of class fitNewDocuments
#' 
#' \item{theta}{a matrix with one row per document contain the document-topic proportions at the posterior mode}
#' \item{eta}{the mean of the variational posterior, only provided when posterior is requested. 
#' Matrix of same dimension as theat}
#' \item{nu}{a list with one element per document containing the covariance matrix of the variational posterior.
#' Only provided when posterior is requested.}
#' \item{phis}{a list with one element per K by V' matrix containing the variational distribution for each token 
#' (where V' is the number of unique words in the given document.  They are in the order of appearence in the document. 
#' For words repeated more than once the sum of the column is the number of times that token appeared.  This is only
#' provided if the posterior is requested.}
#' \item{bound}{a vector with one element per document containing the approximate variational lower bound. This is only
#' provided if the posterior is requested.}
#' \item{beta}{A list where each element contains the unlogged topic-word distribution for each level of the content covariate.
#' This is only provided if prior is requested.}
#' \item{betaindex}{a vector with one element per document indicating which element of the beta list the documents pairs with.
#' This is only provided if prior is requested.}
#' \item{mu}{a matrix where each column includes the K-1 dimension prior mean for each document. This is only provided if prior is requested.}
#' \item{sigma}{a K-1 by K-1 matrix containing the prior covariance. This is only provided if prior is requested.}
#' @seealso \code{\link{alignCorpus}} \code{\link{optimizeDocument}} \code{\link{make.heldout}} \code{\link{makeDesignMatrix}}
#' @examples 
#' #An example using the Gadarian data.  From Raw text to fitted model.
#' #(for a case where documents are all not processed at once see the help
#' # file for alignCorpus)
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' set.seed(02138)
#' #Maximum EM its is set low to make this run fast, run models to convergence!
#' mod.out <- stm(out$documents, out$vocab, 3, prevalence=~treatment + s(pid_rep), 
#'               data=out$meta, max.em.its=5)
#' fitNewDocuments(model=mod.out, documents=out$documents[1:5], newData=out$meta[1:5,],
#'                origData=out$meta, prevalence=~treatment + s(pid_rep),
#'                prevalencePrior="Covariate")
#' @export
fitNewDocuments <- function(model=NULL, documents=NULL, newData=NULL, 
                            origData=NULL, prevalence=NULL, betaIndex=NULL,
                            prevalencePrior=c("Average","Covariate","None"),
                            contentPrior=c("Average","Covariate"),
                            returnPosterior=FALSE,
                            returnPriors=FALSE,
                            designMatrix=NULL, test=TRUE, 
                            verbose=TRUE) {
  
  prevalencePrior <- match.arg(prevalencePrior)
  contentPrior <- match.arg(contentPrior)
  
  #Some checks from the stm() function
  if(!is.list(documents)) stop("documents must be a list, see documentation.")
  if(!all(unlist(lapply(documents, is.matrix)))) stop("Each list element in documents must be a matrix. See documentation.")
  if(any(unlist(lapply(documents, function(x) anyDuplicated(x[1,]))))) {
    stop("Duplicate term indices within a document.  See documentation for proper format.")
  }
  
  #Some additional checks here:
  if(!is.null(newData)) {
    if(nrow(newData)!=length(documents)) stop("Rows in the dataset does not match the length of the new documents")
    if(max(unlist(lapply(documents, function(x) x[1,]))) > length(model$vocab)) stop("Vocabulary in new documents is larger than the original fitted model.")
  }
  
  #set up the counter
  if(verbose) ctevery <- ifelse(length(documents)>100, floor(length(documents)/100), 1)
  
  #Check arguments for prevalencePrior and content prior
   if(!is.null(prevalencePrior)) {
    prevtype <- match.arg(prevalencePrior)
    if(prevtype == "Covariate" & is.null(newData)) stop("Must specify data to use covariate type prevalencePrior")
  } else {
    #if the data is supplied, use it.
    prevtype <- ifelse(is.null(newData), "Average", "Covariate")
  }
  if(!is.null(contentPrior)){
    contenttype <- match.arg(contentPrior)
    if(contenttype == "Covariate" & is.null(betaIndex)) stop("Must specify betaIndex to use covariate type contentPrior")
  } else {
    #if the data is supplied use it.
    contenttype <- ifelse(is.null(betaIndex), "Average", "Covariate") 
  }

  if(is.null(model)) stop("Please specify a fitted STM object for model")
  if(is.null(documents)) stop("Please specify an stm formatted document list for documents argument")
  
  
  K <- model$settings$dim$K
  #Start by Generating mu and sigma
  
  #if None- set to 0's with a broad prior.
  if(prevtype == "None"){
    mu <- matrix(0, nrow=(K-1), ncol=length(documents))
    sigma <- diag(1000, nrow=(K-1))
  }
  #if Average set to the average
  if(prevtype == "Average") {
    #if no covariates just use the actual model parameters
    if(is.null(model$mu$gamma)) {
      mu <- c(model$mu$mu)
      sigma <- model$sigma
    } else {
      #otherwise recalculate the covariance
      mu <- model$mu$mu
      if(ncol(mu)==1) {
        covariance <- crossprod(sweep(model$eta, 2, STATS=as.numeric(mu), FUN="-"))
      } else {
        covariance <- crossprod(model$eta - t(mu)) 
      }
      mu <- rowMeans(mu) #replace with the new averaged mu
      newcovariance <- crossprod(sweep(model$eta, 2, STATS=mu))
      #now to get sigma we just subtract off and replace
      sigma <- model$sigma - covariance/nrow(model$eta) + newcovariance/nrow(model$eta)
    }
  }
  if(prevtype=="Covariate"){
    #Covariates will require a mu and sigma
    if(!is.null(designMatrix)) {
      #if we have manually been passed the design matrix use it here
      if(ncol(model$settings$covariates$X) != ncol(designMatrix)) stop("designMatrix columns do not match original design matrix.")
      X <- designMatrix
    } else {
      #First we need to find the formula
      if(!is.null(prevalence)) {
        #easiest case is that the user provides it.
        if(!inherits(prevalence, "formula")) stop("prevalence argument must provide a formula if specified")
        formula <- prevalence
      } else {
        #okay the user hasn't provided it, let's find it.
        if(!is.null(model$setting$covariates$formula)) {
          #if its a new stm object we've saved it.
          formula <- model$settings$covariates$formula
        } else {
          #third case it isn't in the object, let's reconstruct it from the call
          call <- model$settings$call
          #grab the formula 
          #(code here matches the element of the call, grabs that component, drops the outside parentheses,
          # converts to a formula)
          formula <- as.formula(call[c(match(c("prevalence"), names(call),0L))][[1L]])
        }
      }
     #okay now we have the formula we can call makeDesignMatrix
     X <- makeDesignMatrix(formula, origData, newData, test=test)
    }
    sigma <- model$sigma
    mu <- t(X%*%model$mu$gamma)
  }
  
  #Generate the Content Prior beta and betaindex
  if(contentPrior == "Average") {
    if(length(model$beta$logbeta) == 1) {
      #in this scenario we have the customary case of no content covariates
      betaindex <- rep(1,length(documents))
      beta <- list(exp(model$beta$logbeta[[1]])) #put it in a list we can treat it the same
    } else {
      #in this scenario we have content covariates but no actual data so
      #instead we marginalize.
      betalist <- lapply(model$beta$logbeta, exp)
      weights <- tabulate(model$settings$covariates$betaindex)
      weights <- weights/sum(weights)
      beta <- betalist[[1]]*weights[1]
      for(i in 2:length(betalist)) {
        beta <- beta + betalist[[i]]*weights[i]
      }
      beta <- list(beta) #put it in a list so we can treat them all the same
      #we have now essentially decided to use a single beta
      betaindex <- rep(1,length(documents))
    }
  } else {
    #in this scenario we are going to have different betas
    #based on what the user passed us.
    levels <- model$settings$covariates$yvarlevels
    betaindex <- as.numeric(factor(betaIndex, levels = levels))
    if(any(is.na(betaindex))) stop("NA in the betaIndex variable.  This could be due to the factor names not matching the original model.")
    beta <- lapply(model$beta$logbeta, exp)
  }
  #####
  #Now actually do the optimization
  #####
  #precompute terms for the loop over documents
  sigobj <- try(chol.default(sigma), silent=TRUE)
  if(class(sigobj)=="try-error") {
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
    sigmaentropy <- sum(log(diag(sigobj)))
    siginv <- chol2inv(sigobj)
  }
  
  lambda <- vector(mode="list", length=length(documents))
  if(returnPosterior) {
    nu <- vector(mode="list", length=length(documents))
    phi <- vector(mode="list", length=length(documents))
    bound <- vector(length=length(documents))
  }
  for(i in 1:length(documents)) {
    #update components
    doc <- documents[[i]]
    aspect <- betaindex[i]
    init <- rep(0, K-1)
    #deal with the special case of mu as a vector.
    if(class(mu)!="numeric") {
      mu.i <- mu[,i]
    } else {
      mu.i <- mu 
    }
    #infer the document
    results <- optimizeDocument(doc, eta=init, mu=mu.i, beta=beta[[aspect]],  
                                     sigmainv=siginv, sigmaentropy=sigmaentropy,
                                     posterior=returnPosterior)
    lambda[[i]] <- c(results$lambda)
    if(returnPosterior) {
      nu[[i]] <- results$nu
      phi[[i]] <- results$phi
      bound[i] <- results$bound
    }
    if(verbose && i%%ctevery==0) cat(".")
  }
  lambda <- do.call(rbind, lambda)
  lambda <- cbind(lambda,0)
  theta <- exp(lambda - row.lse(lambda))
  
  #Return the results.  Only use theta if they haven't
  #asked for the priors or posterior.
  toReturn <- list(theta=theta, eta=list(), 
                   nu=list(), phi=list(), 
                   bound=c(), beta=list(),
                   betaindex=c(),
                   mu=c(), sigma=c()) 
  if(returnPosterior) {
    toReturn$nu <- nu
    toReturn$eta <- lambda[,-ncol(lambda), drop=FALSE]
    toReturn$phi <- phi
    toReturn$bound <- bound
  }
    
  if(returnPriors) {
    toReturn$beta <- beta
    toReturn$mu <- mu
    toReturn$sigma <- sigma
    toReturn$betaindex <- betaindex
  }
  class(toReturn) <- "fitNewDocuments"
  return(toReturn)
}
