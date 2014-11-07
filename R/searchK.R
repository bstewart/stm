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
    g$exclus[i]<-mean(unlist(exclusivity(model, M=M, frexw=.7)))
    g$semcoh[i]<-mean(unlist(semanticCoherence(model, heldout$documents, M)))
    g$heldout[i]<-eval.heldout(model, heldout$missing)$expected.heldout    
    g$residual[i]<-checkResiduals(model,heldout$documents)$dispersion
    g$bound[i]<-max(model$convergence$bound)
    g$lbound[i]<-max(model$convergence$bound) + lfactorial(model$settings$dim$K)
    g$em.its[i]<-length(model$convergence$bound)    
  }
  g <- as.data.frame(g)
  toreturn <- list(results=g, call=match.call(expand.dots=TRUE))
  class(toreturn)<- "searchK"
  return(toreturn)
}