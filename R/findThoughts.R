findThoughts <- function(model, texts=NULL, topics=NULL, n=3){
  
  theta <- model$theta
  if(is.null(topics)) topics <- 1:ncol(theta)
  if(!is.null(texts) && length(texts)!=nrow(theta)) stop("Number of provided texts and number of documents modeled do not match")
  x <- apply(theta,2,order,decreasing=TRUE)
  index <- x[1:n,topics, drop=FALSE]
  colnames(index) <- paste("Topic", topics)
  if(is.null(texts)) return(index)
  
  if(is.factor(texts)){
    warning("texts are of type 'factor.'  Converting to character vectors.  Use 'as.character' to avoid this warning in the future.")
    texts <- as.character(texts)  
  }
  
  toprint <- c()
  output <- list()
  j<-topics
  for(i in 1:length(topics)) {
    docs <- paste(c(texts[index[,i]]), collapse="\n \t")
    toprint[i] <- sprintf("\n Topic %i: \n \t %s", j[i], docs)
    output[[i]] <- as.matrix(texts[index[,i]])
    colnames(output[[i]]) <- paste("Topic", i)
  }
  cat(toprint)
  return(list(index=index, docs=output))
}
