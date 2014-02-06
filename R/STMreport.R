
report <- function(model) {
  verbose <- model$settings$verbose
  itnum <- model$convergence$its-1
  vocab <- model$vocab
  if(verbose) {
    convstring <- sprintf("Completing Iteration %i (Approx. Bound = %.3e) \n", itnum, model$convergence$bound[itnum])
    cat(convstring)
  }
  printState <- verbose && itnum%%model$settings$topicreportevery==0
  if(printState) {
    if(!is.null(model$vocab) & model$settings$kappa$LDAbeta) {
      wordmat <- apply(model$beta$logbeta[[1]],1,function(x) vocab[order(x,decreasing=TRUE)[1:5]])
      labs <- apply(wordmat, 2, function(x) paste(x,collapse=", "))
      toprint <- sprintf("Topic %i: %s \n", 1:length(labs), labs)
      cat(toprint)
    }
    if(!is.null(model$vocab) & !model$settings$kappa$LDAbeta) {
      out <- lapply(model$beta$kappa$params, function(x) {
        windex <- order(x,decreasing=TRUE)[1:5]
        windex <- windex[x[windex]>0]
        vocab[windex]
        }) 
      labs <- unlist(lapply(out, function(x) paste(x, collapse=", ")))
      K <- model$settings$dim$K
      A <- model$settings$dim$A
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


