#Make a word cloud.
cloud <- function(stmobj, topic=NULL, type=c("model", "documents"), documents, thresh=.9, max.words=100, ...) {
  if(!requireNamespace("wordcloud", quietly=TRUE)) {
    stop("wordcloud package required to use this function.")
  } else {
    if(class(stmobj)!="STM") stop("cloud function only works for STM models.  See wordcloud package for general tools.")
    if(length(topic)>1) stop("Please only select 1 topic.")
    mod <- stmobj
    type <- match.arg(type)
    vocab <- mod$vocab
    #if they didn't specify a topic overwrite the choice to documents
    if(is.null(topic)) type <- "documents"
    
    if(type=="model") {
      #Here we are interested in the model parameters
      if(length(mod$beta$logbeta)==1) {
        #in the case with no content covariates its simply the reweighted p(w|z)
        vec <- exp(mod$beta$logbeta[[1]])[topic,]*sum(mod$settings$dim$wcounts$x)  
      } else {
        #in the case with content covariates we need to reweight
        levels <- table(mod$settings$covariates$betaindex)
        weights <- levels/sum(levels)
        #now average over the weights
        vec <- weights[1]*exp(mod$beta$logbeta[[1]])[topic,]
        for(i in 2:length(mod$beta$logbeta)) {
          vec <- vec + weights[i]*exp(mod$beta$logbeta[[i]])[topic,]
        }
        #and finally reweight by the marginal word distribution
        vec <- vec*sum(mod$settings$dim$wcounts$x)
      } 
    } else {
      #Here we care about documents with high theta loadings on topic
      if(is.null(topic)) {
        #if no topic is specified its just working with the marginals
        vec <- mod$settings$dim$wcounts$x
      } else {
        if(is.null(documents)) stop("documents needed to give topic specific document values.")
        #Subset to documents fulfulling threshold
        docnums<- which(mod$theta[,topic]>thresh)
        if(length(docnums)==0) stop(sprintf("No documents have a topic loading higher than %s", thresh))
        subdoc <- documents[docnums]
        #Aggregate over the indices to get the margins
        indices <- unlist(lapply(subdoc, "[", 1, ))
        counts <- unlist(lapply(subdoc, "[", 2, ))
        out <- aggregate(counts, by=list(indices), FUN=sum)
        #Fill in the vector of zeroes (to account for words that don't show up in the subpopulation)
        vec <- rep(0, length(vocab))
        vec[out$Group.1] <- out$x
      }
    }
    wordcloud::wordcloud(words=vocab,freq=vec, max.words=max.words, ...)
  }
}