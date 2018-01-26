#' Plots semantic coherence and exclusivity for each topic.
#' 
#' Plots semantic coherence and exclusivity for each topic.  Does not support
#' models with content covariates. 
#' 
#' Each model has semantic coherence and exclusivity values associated with
#' each topic.  This function plots these values and labels each with its topic
#' number.
#' 
#' @param model Output from stm, or a selected model from selectModel.
#' @param documents The documents (see \code{\link{stm}} for format).
#' @param labels Vector of number corresponding to topic numbers.
#' @param M Number of words to use in semantic coherence and exclusivity
#' calculations
#' @param xlab Character string that is x axis title. This should be semantic
#' coherence.
#' @param ylab Character string that is y axis title. This should be
#' exclusivity.
#' @param ...  Other plotting parameters from igraph.
#' @examples
#' 
#' \dontrun{
#' 
#'   #Semantic Coherence calculations require the original documents so we need
#'   #to reconstruct them here.
#'   temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#'   meta<-temp$meta
#'   vocab<-temp$vocab
#'   docs<-temp$documents
#'   out <- prepDocuments(docs, vocab, meta)
#'   docs<-out$documents
#'   vocab<-out$vocab
#'   meta <-out$meta
#'   topicQuality(model=gadarianFit, documents=docs)
#' }
#' @export
topicQuality <- function(model, documents, xlab="Semantic Coherence", ylab="Exclusivity", labels=1:ncol(model$theta), M=10,...){
  # Convert the corpus to the internal STM format
  args <- asSTMCorpus(documents)
  documents <- args$documents

  if(length(model$beta$logbeta)<2){
    semcoh <- semanticCoherence(model,documents=documents, M=M)
    exclusivity <- exclusivity(model, M=M)
    print(semcoh)
    print(exclusivity)
    plot(0, 0, xlab=xlab, ylab=ylab, col="white", xlim=c(min(semcoh), max(semcoh)),
         ylim=c(min(exclusivity), max(exclusivity)),...)
    for(i in 1:length(labels)){
      text(semcoh[i], exclusivity[i], paste("Topic", labels[i]))
    }
  }
  if(length(model$beta$logbeta)>1){
    semcoh <- semanticCoherence(model,documents=documents, M=M)
    for(i in 1:length(labels)){
      print(paste("Topic", i, "has", semcoh[i], "semantic coherence"))
    }
  }
}
