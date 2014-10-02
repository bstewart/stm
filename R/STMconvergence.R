# Simplified convergence checks.
###########
# Changes:
# 1) removed word checks
# 2) distinguish between converged and stopping the em algorithm
convergence.check <- function(bound.ss, convergence, settings) {
  
  #unpack the relevant pieces from the settings
  verbose <- settings$verbose
  emtol <- settings$convergence$em.converge.thresh
  maxits <- settings$convergence$max.em.its
  
  #initialize the convergence object if empty
  if(is.null(convergence)) convergence <- list(bound=c(), its=1, converged=FALSE, stopits=FALSE)

  #fill in the current bound
  convergence$bound[convergence$its] <- sum(bound.ss)

  #if not the first iteration
  if(convergence$its > 1) {
    old <- convergence$bound[convergence$its-1]
    new <- convergence$bound[convergence$its]
    convergence.check <- (new-old)/abs(old)
    #if(convergence.check < emtol & convergence.check > 0) {
    if(convergence.check < emtol) {
      convergence$converged <- TRUE
      convergence$stopits <- TRUE
      if(verbose) cat("Model Converged \n")
      return(convergence)
    }
  }
  
  if(convergence$its==maxits) {
    if(verbose) cat("Model Terminated Before Convergence Reached \n")
    convergence$stopits <- TRUE
    return(convergence)
  }
  convergence$its <- convergence$its + 1
  return(convergence)
}
 


