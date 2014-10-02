#New reporting function.  
# major changes include
# 1) more explicit inputs
# 2) reporting bound/token and/or relative change
# 3) assumption that we are verbose
# 4) assumption that there is some vocab.

report <- function(convergence, ntokens, beta, vocab, topicreportevery, verbose) {
  #compute some components
  itnum <- convergence$its-1 #subtract one because we already iterated the count
  tokenll <- convergence$bound[itnum]/ntokens
  
  #if its the first iteraiton, don't try to print the change
  if(itnum==1) {
    msg <- sprintf("Completing Iteration %i (approx. per word bound = %.3f) \n", 
                   itnum, tokenll) 
  } else {
    old <- convergence$bound[itnum-1]
    new <- convergence$bound[itnum]
    change <- (new-old)/abs(old)
    msg <- sprintf("Completing Iteration %i (approx. per word bound = %.3f, relative change = %.3e) \n", 
                   itnum, tokenll, change) 
  }
  cat(msg)
  
  ####
  # print a summary of the topics if desired.
  
  printState <- verbose && itnum%%topicreportevery==0
  if(printState) {
    if(is.null(beta$kappa)) {
      wordmat <- apply(beta$beta[[1]],1,function(x) vocab[order(x,decreasing=TRUE)[1:5]])
      labs <- apply(wordmat, 2, function(x) paste(x,collapse=", "))
      toprint <- sprintf("Topic %i: %s \n", 1:length(labs), labs)
      cat(toprint)
    } else {
      out <- lapply(beta$kappa$params, function(x) {
        windex <- order(x,decreasing=TRUE)[1:5]
        windex <- windex[x[windex]>0]
        vocab[windex]
        }) 
      labs <- unlist(lapply(out, function(x) paste(x, collapse=", ")))
      K <- nrow(beta$beta[[1]])
      A <- length(beta$beta)
      topics <- sprintf("Topic %i: %s \n", 1:K, labs[1:K])
      cat(topics)
      if(A> 1) {
        i1 <- K + 1
        i2 <- K + A
        aspects <- sprintf("Aspect %i: %s \n", 1:A, labs[i1:i2])
        cat(aspects)
      }
    } 
  }
}


