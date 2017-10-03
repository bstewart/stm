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
#' Due to the need to calculate the heldout-likelihood \code{N} documents have
#' \code{proportion} of the documents heldout at random.  This means that even
#' with the default spectral initialization the results can change from run to run.
#' When the number of heldout documents is low or documents are very short, this also
#' means that the results can be quite unstable.  For example: the \code{gadarian} code
#' demonstration below has heldout results based on only 34 documents and approximately
#' 150 tokens total.  Clearly this can lead to quite disparate results across runs.  By 
#' contrast default settings for the \code{poliblog5k} dataset would yield a heldout sample
#' of 500 documents with approximately 50000 tokens for the heldout sample.  We should expect
#' this to be substantially more stable.
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
#' @param cores Number of CPUs to use for parallel computation
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
                    N = floor(.1*length(documents)), proportion = .5, 
                    heldout.seed = NULL, M = 10, cores = 1, ...) {

  #Make a heldout dataset
  heldout <- make.heldout(documents,vocab, N=N, proportion=proportion, 
                          seed=heldout.seed)

  # warnings
  if( "content" %in% names(list(...)) ) {
    warning("Exclusivity calculation only designed for models without content covariates", call.=FALSE)
  }
  
  # single core
  if (cores == 1) {
      g <- list()
      for (i in seq_along(K)) { # loop produces nicer printout than lapply
          g[[i]] <- get_statistics(K[i], heldout=heldout, init.type=init.type,M=M,...)
      }
  # multi core
  } else {
      cat("Using multiple-cores.  Progress will not be shown. \n")
      g <- parallel::mclapply(K, get_statistics, mc.cores = cores, heldout=heldout, init.type=init.type,
                              M=M,...)
  } 

  # output
  g <- do.call('rbind', g)
  g <- as.data.frame(g)
  toreturn <- list(results=g, call=match.call(expand.dots=TRUE))
  class(toreturn)<- "searchK"
  return(toreturn)
}

# compute statistics for a particular number of topics k
get_statistics <- function(k, heldout, init.type, M, ...) { # k = one particular topic number; K = vector of topic numbers
  out <- NULL # output vector
  out[['K']] <- k
  #run stm
  model <- stm(documents=heldout$documents,vocab=heldout$vocab,
               K=k, init.type=init.type, ...)
  #calculate values to return
  if( !"content" %in% names(list(...)) ) {  # only calculate exclusivity for models without content covariates
    out[['exclus']] <- mean(unlist(exclusivity(model, M=M, frexw=.7)))
    out[['semcoh']] <- mean(unlist(semanticCoherence(model, heldout$documents, M)))
  } 
  out[['heldout']] <- eval.heldout(model, heldout$missing)$expected.heldout    
  out[['residual']] <- checkResiduals(model,heldout$documents)$dispersion
  out[['bound']] <- max(model$convergence$bound)
  out[['lbound']] <- max(model$convergence$bound) + lfactorial(model$settings$dim$K)
  out[['em.its']] <- length(model$convergence$bound)    
  return(out)
}