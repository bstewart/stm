#' Find Thoughts
#' 
#' Outputs most representative documents for a particular topic. Use this in
#' order to get a better sense of the content of actual documents with a high
#' topical content.
#' 
#' Returns the top \code{n} documents ranked by the MAP estimate of the topic's
#' theta value (which captures the modal estimate of the proportion of word
#' tokens assigned to the topic under the model). Setting the \code{thresh}
#' argument allows the user to specify a minimal value of theta for returned
#' documents. Returns document indices and top thoughts.
#' 
#' The \code{plot.findThoughts} function is a shortcut for the \code{plotQuote}
#' function.
#' 
#' @aliases findThoughts print.findThoughts plot.findThoughts
#' @param model Model object created by \code{stm}.
#' @param texts A character vector where each entry contains the text of a
#' document.  Must be in the same order as the documents object.
#' @param topics The topic number or vector of topic numbers for which you want
#' to find thoughts.  Defaults to all topics.
#' @param n The number of desired documents to be displayed per topic.
#' @param thresh Sets a minimum threshold for the estimated topic proportion
#' for displayed documents.  It defaults to imposing no restrictions.
#' @return A \code{findThoughts} object:
#' \item{index}{List with one entry per
#' topic.  Each entry is a vector of document indices.} 
#' \item{docs}{List with
#' one entry per topic.  Each entry is a character vector of the corresponding
#' texts.}
#' @seealso \code{\link{plotQuote}}
#' @examples
#' \dontrun{
#' 
#' findThoughts(gadarianFit, texts=gadarian$open.ended.response, topics=c(1,2), n=3)
#' 
#' #We can plot findThoughts objects using plot() or plotQuote
#' thought <- findThoughts(gadarianFit, texts=gadarian$open.ended.response, topics=1, n=3)
#' 
#' #plotQuote takes a set of sentences
#' plotQuote(thought$docs[[1]])
#' 
#' #we can use the generic plot as a shorthand which will make one plot per topic
#' plot(thought)
#' 
#' #we can select a subset of examples as well using either approach
#' plot(thought,2:3)
#' plotQuote(thought$docs[[1]][2:3])
#' }
#' @export
findThoughts <- function(model, texts=NULL, topics=NULL, n=3, thresh=0.0){
  #Grab up to n texts which are above the threshold
  theta <- model$theta
  if(is.null(topics)) topics <- 1:ncol(theta)
  if(!is.null(texts) && length(texts)!=nrow(theta)) stop("Number of provided texts and number of documents modeled do not match")
  #since we are allowing n to be infinity have to recode if its too high.
  if(n > nrow(theta)) n <- nrow(theta)
  if(n < 1) stop("Must request at least one returned document.")
  
  out <- list()
  for(i in 1:length(topics)) {
    k <- topics[i]
    #grab the values and the rank
    index <- order(theta[,k], decreasing=TRUE)[1:n]
    val <- sort(theta[,k], decreasing=TRUE)[1:n]
    #subset to those values which meet the threshold
    index <- index[which(val>=thresh)]
    out$index[[i]] <- index
    #grab the associated texts
    if(!is.null(texts)) out$docs[[i]] <- texts[index]
  }
  names(out$index) <- paste("Topic", topics)
  if(!is.null(texts)) names(out$docs) <- paste("Topic", topics)
  
  class(out) <- "findThoughts"
  
  if(is.null(texts)) return(out)
  
  out$docs <- lapply(out$docs, unlist)
  if(is.factor(texts)){
    warning("texts are of type 'factor.'  Converting to character vectors.  Use 'as.character' to avoid this warning in the future.")
    out$docs <- lapply(out$docs, as.character)  
  }

  return(out)
}

#' @export
print.findThoughts <- function(x,...) {
  toprint <- vector(length=length(x$docs))
  for(i in 1:length(x$docs)) {
    docs <- paste(x$docs[[i]], collapse="\n \t")
    toprint[i] <- sprintf("\n %s: \n \t %s", names(x[[1]])[i], docs)
  }
  cat(toprint)
}

#' @export
plot.findThoughts <- function(x, sentences=NULL, ...) {
  for(i in 1:length(x$docs)) {
    if(is.null(sentences)) sentences <- 1:length(x$docs[[i]])
    plotQuote(x$docs[[i]][sentences],...)
  }
}