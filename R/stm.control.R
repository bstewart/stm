#Workhorse Function for the STM model
#compared to the original we have more initializations, 
# more explicit options, trimmed fat, memoization

stm.control <- function(documents, vocab, settings, model) {
  
  globaltime <- proc.time()
  verbose <- settings$verbose
  ##########
  #Step 1: Initialize Parameters
  ##########
  ngroups <- settings$ngroups
  if(is.null(model)) {
    if(verbose) cat("Beginning Initialization.\n")
    #initialize
    model <- stm.init(documents, settings)
    #if we were using the Lee and Mimno method of setting K, update the settings
    if(settings$dim$K==0) settings$dim$K <- nrow(model$beta[[1]])
    #unpack
    mu <- list(mu=model$mu)
    sigma <- model$sigma
    beta <- list(beta=model$beta)
    if(!is.null(model$kappa)) beta$kappa <- model$kappa
    lambda <- model$lambda
    convergence <- NULL 
    #discard the old object
    rm(model)
  } else {
    if(verbose) cat("Restarting Model...\n")
    #extract from a standard STM object so we can simply continue.
    mu <- model$mu
    beta <- list(beta=lapply(model$beta$logbeta, exp))
    if(!is.null(model$beta$kappa)) beta$kappa <- model$beta$kappa
    sigma <- model$sigma
    lambda <- model$eta
    convergence <- model$convergence
    #manually declare the model not converged or it will stop after the first iteration
    convergence$stopits <- FALSE
    convergence$converged <- FALSE
    #iterate by 1 as that would have happened otherwise
    convergence$its <- convergence$its + 1 
  }    
  
  #Pull out some book keeping elements
  ntokens <- sum(settings$dim$wcounts$x)
  betaindex <- settings$covariates$betaindex
  stopits <- FALSE
  if(ngroups!=1) {
    groups <- cut(1:length(documents), breaks=ngroups, labels=FALSE) 
  }
  suffstats <- vector(mode="list", length=ngroups)
  
  ############
  #Step 2: Run EM
  ############
  while(!stopits) {
    
    #one set of updates with groups, another without.
    if(ngroups!=1) {
      #####
      #Blocked Updates
      #####
      for(i in 1:ngroups) {
        t1 <- proc.time()
        #update the group id
        gindex <- which(groups==i)
        #construct the group specific sets
        gdocs <- documents[gindex]
        if(is.null(mu$gamma)) {
          gmu <- mu$mu
        } else {
          gmu <- mu$mu[,gindex]
        }
        gbetaindex <- betaindex[gindex]
        glambda <- lambda[gindex,]
        
        #run the model
        suffstats[[i]] <- estep(documents=gdocs, beta.index=gbetaindex, 
                                update.mu=(!is.null(mu$gamma)),  
                                beta$beta, glambda, gmu, sigma, 
                                verbose)
        if(verbose) {
          msg <- sprintf("Completed Group %i E-Step (%d seconds). \n", i, floor((proc.time()-t1)[3]))
          cat(msg)
        }
        t1 <- proc.time()
        
        #if all slots are full.  Combine and run M-step
        if(!any(unlist(lapply(suffstats, is.null)))) {
          #Combine the sufficient statistics
          #(note this is somewhat kludgier than I would prefer
          # but it isn't very costly in terms of time so its fine)
          sigma.ss <- suffstats[[1]]$sigma
          lambda <- suffstats[[1]]$lambda
          beta.ss <- suffstats[[1]]$beta
          bound.ss <- suffstats[[1]]$bound
          for(j in 2:ngroups) {
            sigma.ss <- sigma.ss + suffstats[[j]]$sigma
            lambda <- rbind(lambda, suffstats[[j]]$lambda)
            for(a in 1:length(beta.ss)) {
              beta.ss[[a]] <- beta.ss[[a]] + suffstats[[j]]$beta[[a]]
            }
            bound.ss <- c(bound.ss, suffstats[[j]]$bound)
          }
          # Now do the updates themselves
          mu <- opt.mu(lambda=lambda, mode=settings$gamma$mode, 
                       covar=settings$covariates$X, settings$gamma$enet)
          sigma <- opt.sigma(nu=sigma.ss, lambda=lambda, 
                             mu=mu$mu, sigprior=settings$sigma$prior)
          beta <- opt.beta(beta.ss, beta$kappa, settings)
 
          if(verbose) {
           #M-step message
            timer <- floor((proc.time()-t1)[3])
            msg <- ifelse(timer>1, 
                          sprintf("Completed M-Step (%d seconds). \n", floor((proc.time()-t1)[3])),
                          "Completed M-Step. \n") 
            cat(msg)
          }
        }
      }
    } else {
      #####
      # Non-Blocked Updates
      #####
      t1 <- proc.time()
      #run the model
      suffstats <- estep(documents=documents, beta.index=betaindex, 
                              update.mu=(!is.null(mu$gamma)),  
                              beta$beta, lambda, mu$mu, sigma, 
                              verbose)
      msg <- sprintf("Completed E-Step (%d seconds). \n", floor((proc.time()-t1)[3]))
      if(verbose) cat(msg)
      t1 <- proc.time()
      sigma.ss <- suffstats$sigma
      lambda <- suffstats$lambda
      beta.ss <- suffstats$beta
      bound.ss <- suffstats$bound
      #do the m-step
      mu <- opt.mu(lambda=lambda, mode=settings$gamma$mode, 
                   covar=settings$covariates$X, settings$gamma$enet)
      sigma <- opt.sigma(nu=sigma.ss, lambda=lambda, 
                         mu=mu$mu, sigprior=settings$sigma$prior)
      beta <- opt.beta(beta.ss, beta$kappa, settings)
      if(verbose) {
        timer <- floor((proc.time()-t1)[3])
        msg <- ifelse(timer>1, 
                      sprintf("Completed M-Step (%d seconds). \n", floor((proc.time()-t1)[3])),
                      "Completed M-Step. \n") 
        cat(msg)
      }
    }
    #Convergence
    convergence <- convergence.check(bound.ss, convergence, settings)
    stopits <- convergence$stopits

    #Print Updates if we haven't yet converged
    if(!stopits & verbose) report(convergence, ntokens=ntokens, beta, vocab, 
                                       settings$topicreportevery, verbose)
  }
  #######
  #Step 3: Construct Output
  #######
  time <- (proc.time() - globaltime)[3]
  #convert the beta back to log-space
  beta$logbeta <- beta$beta
  for(i in 1:length(beta$logbeta)) {
    beta$logbeta[[i]] <- log(beta$logbeta[[i]])
  }
  beta$beta <- NULL
  lambda <- cbind(lambda,0)
  model <- list(mu=mu, sigma=sigma, beta=beta, settings=settings,
                vocab=vocab, convergence=convergence, 
                theta=exp(lambda - row.lse(lambda)), 
                eta=lambda[,-ncol(lambda), drop=FALSE],
                invsigma=solve(sigma), time=time, version=utils::packageDescription("stm")$Version)
  class(model) <- "STM"  
  return(model)
}



    


