#TODO: Create an exported version of the design matrix function
#      Update parseFormulas to use the new thing
#      Write a testing function
#      Content covariates
#      Write the part that actually does the fitting
#      Write the return structure and document
#      Write an example
#      Put in seealso links



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
#' \strong{Testing:} The complication of this function is getting the transformations right.  Functions
#' like \code{s()}, \code{bs()}, \code{ns()} and \code{poly()} produce an output that is dependent on
#' all the observations passed to it.  Thus when doing the fitting we need to ensure the new prediction
#' matrix is created in the same way as when the model was original fit.  \code{fitNewDocuments} will
#' attempt to solve this problem internally but we can't solve every edge case. As such we have built-in
#' a way of testing if the process had worked.  As long as \code{test=TRUE} the function will take subsets
#' of the original data and evaluate if it has succeeded in reconstructing those values.  If the values
#' produced are different it will issue an error and stop.
#' 
#' \strong{Passing a Design Matrix}  Advanced users may wish to circumvent this process and pass their
#' own design matrix possibly because they used their own function for transforming the original input
#' variables.  This can be done by passing the design matrix using \code{designMatrix}.  The columns
#' need to match the ordering of the design matrix for the original \code{stm} object.  The design matrix
#' in an stm model called \code{stmobj} can be found in \code{stmobj$settings$covariates$X} which can
#' in turn be used to check that you have formatted your result correctly. If you want even more fine-grained
#' control we recommend you directly use the optimization function \code{\link{optimizeDocument}}
#' 
#' \strong{specialFunctions} As mentioned above the key difficulty is functions in the formulas which
#' are dependent for their form on the data passed to them (like splines e.g \code{bs()}, but unlike 
#' e.g. \code{log()}).  This function allows for the special functions \code{"s", "ns", "bs", "poly", "pspline"}
#' internally what we do is create a temporary function that calls the relevant predict function
#' e.g. \code{predict.bs()} using the generic assigned to that function and the original data.  If you
#' want to use a different function which also has a predict method associated with it, you can pass it
#' here and it will be incorporated.
#' 
#' @param model the originally fit STM object.
#' @param documents the new documents to be fit. These documents must be in the stm format and
#'  be numbered in the same way as the documents in the original model with the same dimension of vocabulary.
#' @param newData the metadata for the prevalence prior which goes with the unseen documents. As in
#' the original data this cannot have any missing data.
#' @param origData the original metadata used to fit the STM object. 
#' @param prevalence the original formula passed to prevalence when \code{stm()} was called. The function
#' will try to reconstruct this.
#' @param betaIndex an integer vector which indicates which level of the content covariate is used
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
#' @param test a test of the functions ability to reconstruct the original functions.
#' @param designMatrix an option for advanced users to pass an already constructed design matrix for
#' prevalence covariates.  This will override the options in \code{newData}, \code{origData} and
#' \code{test}.  See details below- please do not attempt to use without reading carefully.
#' @param specialFunctions an option for advanced users.  This allows you to pass a character
#' vector containing the names of any special functions that you may wish to use on the data
#' inside your formula. See details below.  By default we include \code{"s", "ns", "bs", "poly", "pspline"}
#' @examples
#' 
#' @seealso \code{\link{optimizeDocument}} \code{\link{make.heldout}}
#' @export
fitNewDocuments <- function(model, documents, newData=NULL, 
                            origData=NULL, prevalence=NULL, betaIndex=NULL,
                            prevalencePrior=c("Average","Covariate","None"),
                            contentPrior=c("Average","Covariate"),
                            returnPosterior=FALSE,
                            returnPriors=FALSE,
                            test=TRUE, designMatrix=NULL,
                            specialFunctions=NULL) {
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
  if(prevtype=="Covariate"){
    #Covariates will require a mu and sigma
    if(!missing(designMatrix)) {
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
      
      origData$lpid_rep <- asinh(origData$pid_rep)
      formula <- ~treatment*s(pid_rep) + s(lpid_rep)
      
      #Second, we now have the formula, let's evaluate in the context of the original data
      termobj <- terms(formula, data=origData)
      mf <- model.frame(termobj, data=origData)
    
      #now let's identify the specials and build out functions for them
      elements <- colnames(attr(terms(formula), "factors"))
      specialtypes <- c("s", "ns", "bs", "poly", "pspline")
      if(!is.null(specialFunctions)) specialtypes <- c(specialtypes,specialFunctions)
      special <- attr(terms(formula, 
                            specials=specialtypes),
                      "special")
      newformula <- formula
      for(i in 1:length(special)) {
        #for each special check if there are any
        if(is.null(special[[i]])) next
        #okay now that we know there is at least one
        #now we step through five steps
        
        #1) loop over elements of specials
        for(j in 1:length(special[[i]])) {
          #2) find the element of interest
          elem <- elements[special[[i]][j]]
          #3) create a function name based on names(special)[i].  
          funcname <- sprintf("%sfn%i",names(special)[i],j)
          #4) write a predict function using that name
          #   this calls predict on the newdata (x) using the old model frame
          assign(funcname, function(x,...) predict(mf[[elem]],x))
          #5) figure out what just the variable name is.  There is probably a more elegant way to do this
          # the more obvious way is via all.vars() but this doesn't line us up with elements
          # also get_all_vars but same problem.
          # The basic strategy is to use reverse to get only the first and last parenthesis so functions can be nested
          varname <- stringi::stri_split_fixed(elem, "(", n=2)[[1]][2]
          varname <- stringi::stri_reverse(varname)
          varname <- stringi::stri_split_fixed(varname, ")",n=2)[[1]][2]
          varname <- stringi::stri_reverse(varname)
          #6) replace the full element in the formula with the newfunction called on the original variable
          replace <- sprintf("%s(%s)", funcname, varname)
          newformula <- as.character(newformula)
          newformula[2] <- stringi::stri_replace_all_fixed(as.character(newformula)[2],
                                                           replacement=replace,
                                                           pattern=elem)
          newformula <- as.formula(newformula) #coerce back to formula
          #7) move on to the next one
        }
      }
      #now we have a new variable called newformula which has an updated formula for us to call!   
      X <- try(Matrix::sparse.model.matrix(termobj,data=origData),silent=TRUE)
      if(class(X)=="try-error") stop("Error creating model matrix. Check your formula.")
      
      if(test) {
        #write a test here. First I need to export.
      }
    }
    sigma <- model$sigma
    mu <- t(X%*%gamma)
  }
  
  #Generate the Content Prior Matrix
     #the trick here is using the yvarlevels to work back what the factor is.  
     #Should have three types: integer, factor, character
  
  
  
  #Loop Over and fit them all
  
  #Return some structure including the fits.
   
}
