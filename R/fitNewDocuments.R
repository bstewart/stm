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
#' set.  This does not require the unseen documents to observe the covariates.  In the case of 
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
#' @param data the metadata which goes with the unseen documents.
#' @param prevalencePrior three options described in detail below.
#' @param contentPrior two options described in detail below.
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
fitNewDocuments <- function(model, documents, data=NULL,
                            prevalencePrior=c("Average","Covariate","None"),
                            contentPrior=c("Average","Covariate"),
                            returnPosterior=FALSE,
                            returnPriors=FALSE) {
  
  #TODO: argument matching for prevalence and content Priors
  #      decide on defaults and be clear about them
  
  if(missing(model)) stop("Please specify a fitted STM object for model")
  if(missing(documents)) stop("Please specify an stm formatted document list for documents argument")

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
  
 #TODO: write some checks for nonsensical combinations like covariates without specifying data or
  #     contentPrior models
  
  
 #Generate the Prevalence Prior Matrix
    #if None- set to 0's with a broad prior.

   #Generate the Content Prior Matrix
    #
  
  #Loop Over and fit them all
  
  #Return some structure including the fits.
   
}