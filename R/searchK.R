searchK <- function(documents, vocab, K, init.type = "Spectral", 
                    N=floor(.1*length(documents)), proportion=.5, heldout.seed=NULL,
                    M=10,...) {

  gnames<-c("K","heldout","residual","bound","lbound","exclus","semcoh","em.its")
  g<-vector(mode="list", length=length(gnames))
  names(g) <- gnames
  for(i in 1:length(g)) {
    g[[i]] <- vector(length=length(K))
  }
  
  heldout <- make.heldout(documents,vocab, N=N, proportion=proportion, seed=heldout.seed)
  docs <- heldout$documents
  vocab <- heldout$vocab
  
  for(i in 1:length(K)) {
    
    g$K[i]<-K[i]
    
        if(init.type!="Spectral"){
        models<- selectModel(documents=documents, vocab=vocab, K=K[i], init.type=init.type, ...)
        j<-paretosingle(models)
        storage<-models$runout[[j]]
        g$exclus[i]<-storage$exclusivity
        g$semcoh[i]<-storage$semcoh
        j<-NULL
        
    
        g$heldout[i]<-eval.heldout(storage, heldout$missing)$expected.heldout
        g$residual[i]<-checkResiduals(storage,docs)$dispersion
        g$bound[i]<-max(storage$convergence$bound)
        g$lbound[i]<-max(storage$convergence$bound) + lfactorial(storage$settings$dim$K)
        
      }
  
        if(init.type=="Spectral"){
          storage<- stm(documents=documents, vocab=vocab, K=K[i], init.type=init.type,...)
          g$exclus[i]<-mean(unlist(exclusivity(storage, M=M, frexw=.7)))
          g$semcoh[i]<-mean(unlist(semanticCoherence(storage, documents, M)))
          g$heldout[i]<- eval.heldout(storage, heldout$missing)$expected.heldout
          if(!is.finite(g$heldout[i])) browser()
          g$residual[i]<-checkResiduals(storage,docs)$dispersion
          g$bound[i]<-max(storage$convergence$bound)  
          g$lbound[i]<-lfactorial(storage$settings$dim$K)+max(storage$convergence$bound)  
          g$em.its[i]<-length(storage$converge$bound)
          
        }
  

  }
  rm(storage)
  g <- as.data.frame(g)
  toreturn <- list(results=g, call=match.call(expand.dots=TRUE))
  class(toreturn)<-c("searchK")
  return(toreturn)
}