#function to choose the pareto dominant model, or randomly choose among model candidates that are not weakly dominanted by other models
#utility function for manyTopics

paretosingle <- function(z) {
  
  m<-matrix(NA,nrow=length(z$semcoh),ncol=2)
  
  for(i in 1:nrow(m)){
    
    m[i,1]<-as.numeric(mean(unlist(z$semcoh[i])))
    
    if(!z$exclusivity[[1]][1]=="Exclusivity not calculated for models with content covariates"){
      m[i,2]<-as.numeric(mean(unlist(z$exclusivity[i])))
    } else {
      stop("manyTopics function not yet designed for models with content variable.") 
    }
  }
 
  d1max <- max(m[,1])
  d2max <- max(m[,2])
  weakcandidates <- m[,1]==d1max | m[,2]==d2max
  strongcandidates <- m[,1]==d1max & m[,2]==d2max
  s = which(strongcandidates)
  w = which(weakcandidates)
  if (length(s)>0) {
    x = s
  }
  else {
    x = w
  }
  if (length(x)==1) {return(x)}
  else {
    return(sample(x,size=1))
  }
}

#' Performs model selection across separate STM's that each assume different
#' numbers of topics.
#' 
#' Works the same as selectModel, except user specifies a range of numbers of
#' topics that they want the model fitted for. For example, models with 5, 10,
#' and 15 topics.  Then, for each number of topics, selectModel is run multiple
#' times. The output is then processed through a function that takes a pareto
#' dominant run of the model in terms of exclusivity and semantic coherence. If
#' multiple runs are candidates (i.e., none weakly dominates the others), a
#' single model run is randomly chosen from the set of undominated runs. 
#' 
#' Does not work with models that have a content variable (at this point).
#' 
#' @param documents The documents to be modeled.  Object must be a list of with
#' each element corresponding to a document.  Each document is represented as
#' an integer matrix with two rows, and columns equal to the number of unique
#' vocabulary words in the document.  The first row contains the 1-indexed
#' vocabulary entry and the second row contains the number of times that term
#' appears.
#' 
#' This is similar to the format in the \pkg{lda} package except that
#' (following R convention) the vocabulary is indexed from one. Corpora can be
#' imported using the reader function and manipulated using the
#' \code{\link{prepDocuments}}.
#' @param vocab Character vector specifying the words in the corpus in the
#' order of the vocab indices in documents. Each term in the vocabulary index
#' must appear at least once in the documents.  See
#' \code{\link{prepDocuments}} for dropping unused items in the vocabulary.
#' @param K A vector of positive integers representing the desired number of
#' topics for separate runs of selectModel.
#' @param prevalence A formula object with no response variable or a matrix
#' containing topic prevalence covariates.  Use \code{s()}, \code{ns()} or
#' \code{bs()} to specify smooth terms. See details for more information.
#' @param content A formula containing a single variable, a factor variable or
#' something which can be coerced to a factor indicating the category of the
#' content variable for each document.
#' @param runs Total number of STM runs used in the cast net stage.
#' Approximately 15 percent of these runs will be used for running a STM until
#' convergence.
#' @param data Dataset which contains prevalence and content covariates.
#' @param init.type The method of initialization.  See \code{\link{stm}}.
#' @param seed Seed for the random number generator. \code{stm} saves the seed
#' it uses on every run so that any result can be exactly reproduced.  When
#' attempting to reproduce a result with that seed, it should be specified
#' here.
#' @param max.em.its The maximum number of EM iterations.  If convergence has
#' not been met at this point, a message will be printed.
#' @param emtol Convergence tolerance.
#' @param verbose A logical flag indicating whether information should be
#' printed to the screen.
#' @param frexw Weight used to calculate exclusivity
#' @param net.max.em.its Maximum EM iterations used when casting the net
#' @param netverbose Whether verbose should be used when calculating net
#' models.
#' @param M Number of words used to calculate semantic coherence and
#' exclusivity.  Defaults to 10.
#' @param \dots Additional options described in details of stm.
#' @return 
#' 
#' \item{out}{List of model outputs the user has to choose from.  Take
#' the same form as the output from a stm model.} 
#' \item{semcoh}{Semantic
#' coherence values for each topic within each model selected for each number
#' of topics.} 
#' \item{exclusivity}{Exclusivity values for each topic within each
#' model selected.  Only calculated for models without a content covariate.}
#' @examples
#' 
#' \dontrun{
#' 
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' meta<-temp$meta
#' vocab<-temp$vocab
#' docs<-temp$documents
#' out <- prepDocuments(docs, vocab, meta)
#' docs<-out$documents
#' vocab<-out$vocab
#' meta <-out$meta
#' 
#' set.seed(02138)
#' storage<-manyTopics(docs,vocab,K=3:4, prevalence=~treatment + s(pid_rep),data=meta, runs=10)
#' #This chooses the output, a single run of STM that was selected,
#' #from the runs of the 3 topic model
#' t<-storage$out[[1]]
#' #This chooses the output, a single run of STM that was selected,
#' #from the runs of the 4 topic model
#' t<-storage$out[[2]]
#' #Please note that the way to extract a result for manyTopics is different from selectModel.
#' }
#' @export
manyTopics <- function(documents, vocab, K, prevalence=NULL, content=NULL, 
                       data = NULL,max.em.its = 100, verbose = TRUE, 
                       init.type = "LDA", 
                       emtol = 1e-05, seed = NULL, runs = 50, 
                       frexw = 0.7, net.max.em.its = 2, 
                       netverbose = FALSE, M = 10,...) {
  out<-list()
  semcoh<-exclusivity<-list()
  for(i in 1:length(K)) {
    
    models<- selectModel(documents, vocab, K[i], prevalence, content, data = data, 
                         max.em.its = max.em.its, verbose = verbose, init.type = init.type, emtol = emtol, seed = seed, runs = runs, 
                         frexw = frexw, net.max.em.its = net.max.em.its, netverbose = netverbose, M = M, 
                         ...)
    j<-paretosingle(models)

    out[[i]]<-models$runout[[j]]
    exclusivity[[i]]<-models$exclusivity[[j]]
    semcoh[[i]]<-models$semcoh[[j]]
    j<-NULL
  }
  
  toreturn<-list(out=out,exclusivity=exclusivity,semcoh=semcoh)
  return(toreturn)
}