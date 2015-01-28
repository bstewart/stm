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

print.findThoughts <- function(x,...) {
  toprint <- vector(length=length(x$docs))
  for(i in 1:length(x$docs)) {
    docs <- paste(x$docs[[i]], collapse="\n \t")
    toprint[i] <- sprintf("\n %s: \n \t %s", names(x[[1]])[i], docs)
  }
  cat(toprint)
}

plot.findThoughts <- function(x, sentences=NULL, ...) {
  for(i in 1:length(x$docs)) {
    if(is.null(sentences)) sentences <- 1:length(x$docs[[i]])
    plotQuote(x$docs[[i]][sentences],...)
  }
}

