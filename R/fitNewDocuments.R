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
#' For each case we offer three choices (and may offer more in the future):
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
#' @param model the originally fit STM object.
#' @param documents these documents must be in the stm format and be numbered in the same way
#' as the documents in the original model with the same dimension of vocabulary!
#' @param data the metadata for the prevalence prior which goes with the unseen documents
#' @param betaIndex an integer vector which indicates which level of the content covariate is used
#' for each unseen document.
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
#' @examples
#' 
#' @export
fitNewDocuments <- function(model, documents, data=NULL, betaIndex=NULL,
                            prevalencePrior=c("Average","Covariate","None"),
                            contentPrior=c("Average","Covariate"),
                            returnPosterior=FALSE,
                            returnPriors=FALSE) {
  #Some checks from the stm() function
  if(!is.list(documents)) stop("documents must be a list, see documentation.")
  if(!all(unlist(lapply(documents, is.matrix)))) stop("Each list element in documents must be a matrix. See documentation.")
  if(any(unlist(lapply(documents, function(x) anyDuplicated(x[1,]))))) {
    stop("Duplicate term indices within a document.  See documentation for proper format.")
  }
  
  #Some additional checks here:
  if(!is.null(data)) {
    if(nrow(data)!=length(documents)) stop("Rows in the dataset does not match the length of the new documents")
    if(max(unlist(lapply(documents, function(x) x[1,]))) > length(model$vocab)) stop("Vocabulary in new documents is larger than the original fitted model.")
  }
  
  
   if(!missing(prevalencePrior)) {
    prevtype <- match.arg(prevalencePrior)
    if(prevtype == "Covariate" & is.null(data)) stop("Must specify data to use covariate type prevalencePrior")
  } else {
    #if the data is supplied, use it.
    prevtype <- ifelse(missing(data), "Average", "Covariate")
  }
  if(!missing(contentPrior)){
    contenttype <- match.arg(contentPrior)
    if(contenttype == "Covariate" & is.null(betaIndex)) stop("Must specify betaIndex to use covariate type contentPrior")
  } else {
    #if the data is supplied use it.
    contenttype <- ifelse(missing(betaIndex), "Average", "Covariate") 
  }

  if(missing(model)) stop("Please specify a fitted STM object for model")
  if(missing(documents)) stop("Please specify an stm formatted document list for documents argument")
  
  
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
        covariance <- crossprod(lambda-t(mu)) 
      }
      mu <- rowMeans(mu) #replace with the new averaged mu
      newcovariance <- crossprod(sweep(model$eta, 2, STATS=mu))
      #now to get sigma we just subtract off and replace
      sigma <- model$sigma - covariance/nrow(model$eta) + newcovariance/nrow(model$eta)
    }
  }
  if(prevtype=="Covariate") {
    #Covariates
    #this is the tough case.  Need to do our best to reconstruct the design matrix for the
    
    #TODO: require old data, modify parseFormulas to allow all basis functions.  Consider
    #including poly in there too.  It has a predict function.  I think ns() and bs() do too
    #so maybe that's the whole gambit?
    
    #I think we should consider exporting the function too...
    str(model)
  }

   #Generate the Content Prior Matrix
  
  #Loop Over and fit them all
  
  #Return some structure including the fits.
   
}