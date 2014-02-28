#Workhorse Function for the Structural Topic Model.
#mostly calls other functions, but this is the high level overview.

stm.control <- function(documents, vocab=NULL, settings) {
  
  #Step 1: Initialize Parameters
  model <- stm.initialize(documents, settings) 
  model$vocab <- vocab #put it in even if null
  history <- list()
  converged <- FALSE
  suffstats <- NULL
  #Step 2: Run EM
  while(!converged) {
    #Save history if we are doing that
    if(settings$keepHistory) {
      history[[(length(history)+1)]] <- list(model=model,suffstats=suffstats)
    }
    t1 <- proc.time()
    #Run the E-Step    
    suffstats <- estep.LN(documents, settings$covariates$betaindex, 
                          logbeta=model$beta$logbeta, 
                          mu=model$mu$mu, sigma=model$sigma, 
                          settings$verbose, lambdacurrent=suffstats$lambda) 
    if(settings$verbose) cat(sprintf("Completed E-Step (%d seconds). \n", floor((proc.time()-t1)[3])))
    t1 <- proc.time()
    #M-Step
    model$mu <- opt.mu(lambda=suffstats$lambda, mode=settings$gamma$mode, 
                       covar=settings$covariates$X, settings$gamma$enet)

    model$sigma <- opt.sigma(nu=suffstats$sigma, lambda=suffstats$lambda, 
                             mu=model$mu$mu, sigprior=settings$sigma$prior)
    
    model$beta <- opt.beta(suffstats$beta, 
                           model$beta$kappa,
                           settings)
    #if verbose- write out something- it just looks silly if it took less than 2 seconds
    #though so don't bother.
    if(settings$verbose) {
      timer <- floor((proc.time()-t1)[3])
      if(timer>1) {
        cat(sprintf("Completed M-Step (%d seconds). \n", floor((proc.time()-t1)[3])))
      } else {
        cat("Completed M-Step. \n")
      }
    }
    
    #Convergence
    model$convergence <- convergence.check(suffstats$bound, model, settings) #assess convergence
    converged <- model$convergence$converged
    
    #Print Updates if we haven't yet converged
    if(!converged) {
      report(model)
    }
  }
  
  #Step 3: Construct Output
  return(make.output(model, suffstats, history))
}

make.output <- function(model,suffstats, history) {
  if(!model$settings$keepHistory) {
    history <- "History Not Kept.  Set keepHistory to true"
  }
  
  eta <- cbind(suffstats$lambda,0)
  model$theta <- exp(eta - row.lse(eta))
  model$eta <- suffstats$lambda
  model$invsigma <- solve(model$sigma)  
  model$history <- history
  class(model) <- "STM"
  return(model)  
}

    


