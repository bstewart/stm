# Functions for checking Convergence in the Structural Topic Model

convergence.check <- function(bound.ss, model, settings) {
  convergence <- model$convergence
  verbose <- settings$verbose
  emtol <- settings$convergence$em.converge.thresh
  maxits <- settings$convergence$max.em.its
  nwords <- min(settings$convergence$topwords,5)
  nwordits <- settings$convergence$topwordits
 
  gettopwords <- function(logbeta, M) {
     newbeta <- do.call(rbind,logbeta)
     apply(newbeta,1,order,decreasing=TRUE)[1:M,]
  }
  
  if(is.null(convergence)) {
    convergence <- list(bound=c(), its=1, wordchangectr=0,converged=FALSE)
  }
  newwords <- gettopwords(model$beta$logbeta, nwords)
  convergence$bound[convergence$its] <- sum(bound.ss)

  if(convergence$its > 1) {
    oldwords <- convergence$priortopwords
    same <- vector(length=ncol(newwords))
    for(i in 1:ncol(newwords)) {
      same[i] <- all(oldwords[,i] %in% newwords[,i])
    }
    if(all(same)) {
      convergence$wordchangectr <- convergence$wordchangectr + 1
    } else {
      convergence$wordchangectr <- 0
    }
    
    old <- convergence$bound[convergence$its-1]
    new <- convergence$bound[convergence$its]
    convergence.check <- (new-old)/abs(old)
    
    if(convergence.check < emtol) {
      convergence$converged <- TRUE
      if(verbose) cat("Model Converged \n")
      return(convergence)
    }
    if(nwordits != 0 & nwordits==convergence$wordchangectr) {
      convergence$converged <- TRUE
      if(verbose) {
        convstring <- sprintf("Model converged (Approx. Bound = %.3e) due to %i consecutive iterations with no change in top %i words within each topic. \n", 
                              new, convergence$wordchangectr, nwords)
        cat(convstring)
      }
      return(convergence)
    }
  }
  
  if(convergence$its==maxits) {
    if(verbose) cat("Model Terminated Before Convergence Reached \n")
    convergence$converged <- TRUE
    return(convergence)
  }
  convergence$priortopwords <- newwords
  convergence$its <- convergence$its + 1
  return(convergence)
}
 


