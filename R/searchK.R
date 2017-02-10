#' Computes diagnostic values for models with different values of K (number of
#' topics).
#' 
#' With user-specified initialization, this function runs selectModel for
#' different user-specified topic numbers and computes diagnostic properties
#' for the returned model. These include exclusivity, semantic coherence,
#' heldout likelihood, bound, lbound, and residual dispersion.
#' 
#' See the vignette for interepretation of each of these measures.  Each of
#' these measures is also available in exported functions:
#' \describe{
#' \item{exclusivity}{\code{\link{exclusivity}}}
#' \item{semantic coherence}{\code{\link{semanticCoherence}}}
#' \item{heldout likelihood}{\code{\link{make.heldout}} and \code{\link{eval.heldout}}}
#' \item{bound}{calculated by \code{\link{stm}} accessible by \code{max(model$convergence$bound)}}
#' \item{lbound}{a correction to the bound that makes the bounds directly comparable \code{max(model$convergence$bound) + lfactorial(model$settings$dim$K)}}
#' \item{residual dispersion}{\code{\link{checkResiduals}}}
#' }
#' 
#' @param documents The documents to be used for the stm model
#' @param vocab The vocabulary to be used for the stmmodel
#' @param K A vector of different topic numbers
#' @param init.type The method of initialization. See \code{\link{stm}} for
#' options.  Note that the default option here is different from the main
#' function.
#' @param N Number of docs to be partially held out
#' @param proportion Proportion of docs to be held out.
#' @param heldout.seed If desired, a seed to use when holding out documents for
#' later heldout likelihood computation
#' @param M M value for exclusivity computation
#' @param ...  Other diagnostics parameters.
#' @return \item{exclus}{Exclusivity of each model.} \item{semcoh}{Semantic
#' coherence of each model.} \item{heldout}{Heldout likelihood for each model.}
#' \item{residual}{Residual for each model.} \item{bound}{Bound for each
#' model.} \item{lbound}{lbound for each model.} \item{em.its}{Total number of
#' EM iterations used in fiting the model.}
#' @seealso \code{\link{plot.searchK}} \code{\link{make.heldout}}
#' @examples
#' 
#' \dontrun{
#' 
#' K<-c(5,10,15) 
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' documents <- out$documents
#' vocab <- out$vocab
#' meta <- out$meta
#' set.seed(02138)
#' K<-c(5,10,15) 
#' kresult <- searchK(documents, vocab, K, prevalence=~treatment + s(pid_rep), data=meta)
#' plot(kresult)
#' 
#' }
#'  
#' @export
searchK <- function(documents, vocab, K, init.type = "Spectral", 
                    N=floor(.1*length(documents)), proportion=.5, 
                    heldout.seed=NULL,
                    M=10,...) {
  #Set up the object to return
  g <- rep(list(vector(length=length(K))), 8)
  names(g) <- c("K","heldout","residual","bound","lbound","exclus","semcoh","em.its")
  
  #Make a heldout dataset
  heldout <- make.heldout(documents,vocab, N=N, proportion=proportion, 
                          seed=heldout.seed)
  
  #Loop over each of the number of topics
  for(i in 1:length(K)) {
    g$K[i]<-K[i]
    #run stm
    model <- stm(documents=heldout$documents,vocab=heldout$vocab,
                 K=K[i], init.type=init.type,...)
    #calculate values to return
    if( !"content" %in% names(list(...)) ) {  # only calculate exclusivity for models without content covariates
      g$exclus[i]<-mean(unlist(exclusivity(model, M=M, frexw=.7)))
      g$semcoh[i]<-mean(unlist(semanticCoherence(model, heldout$documents, M)))
    }
    g$heldout[i]<-eval.heldout(model, heldout$missing)$expected.heldout    
    g$residual[i]<-checkResiduals(model,heldout$documents)$dispersion
    g$bound[i]<-max(model$convergence$bound)
    g$lbound[i]<-max(model$convergence$bound) + lfactorial(model$settings$dim$K)
    g$em.its[i]<-length(model$convergence$bound)    
  }
  g <- as.data.frame(g)
  if( "content" %in% names(list(...)) ) {
    warning("Exclusivity calculation only designed for models without content covariates", call.=FALSE)
    g$exclus <- NULL
    g$semcoh <- NULL
  }
  toreturn <- list(results=g, call=match.call(expand.dots=TRUE))
  class(toreturn)<- "searchK"
  return(toreturn)
}