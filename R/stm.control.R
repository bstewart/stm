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
    # need to test if documents and vocab are the same
    # What to copy:
    #   mu - a topic * document matrix
    #   lambda = eta  - a document * topic matrix
    #   beta -  a data structure relating topics and terms
    if (verbose) cat("Checking if documents or vocabular have changed...\n")
    if (! identical(colnames(mu$mu), names(documents))) {
      cat("Documents not identical\n")
      # 1. Where documents have been removed, delete the entries in mu and lambda/eta
      mu <- model$mu[,-which(! colnames(mu$mu) %in% names(documents) )]
      lambda <- model$eta[-which(! rownames(model$eta) %in% names(documents)),]
      # 2. Reorder documents to match mu and lambda/eta, with new documents at the end
      docorder <- match(names(documents), colnames(mu), nomatch = 0)
      new_docs <- which(docorder == 0)
      new_doc_idxs <- NULL
      if (length(new_docs) > 0) {
        new_doc_idxs <- seq(max(docorder):(max(docorder) + length(new_docs)))
        docorder[new_docs] <- new_doc_idxs
        documents <- documents[docorder]
      }
      # 3. Create new entries in mu and lambda/eta for the new documents
      if (! is.null(new_doc_idxs)) {
        mu <- cbind(mu, matrix(rep(NA,length(new_doc_idxs) * nrow(mu), 
                                   nrow = nrow(mu))))
        lambda <- rbind(lambda, matrix(rep(NA,length(new_doc_idxs) * ncol(lambda), 
                                           ncol = ncol(lambda))))
      }
      # 4. We will have to initialize those entries using new_doc_idxs after lining up and initializing terms
    } else {
      cat("Documents identical\n")
      mu <- model$mu # Rows are topics, columns are documents
      lambda <- model$eta #Rows are documents, columns are topics
    }

    if (! identical(vocab, model$vocab)) {
      cat("Vocabulary not identical\n")
      # beta
      #   beta$kappa - list
      #     beta$kappa$m -- numeric vector for terms
      #     beta$kappa$params -- list of term vectors.  What are the list elements?  Unknown.
      #   beta$nlambda - single number, not copied
      #   beta$logbeta -- list of three matrices
      #     beta$logbeta[1:3] -- each a matrix of K + 1 rows, terms are columns
      # Here don't want to remove entries for removed terms.  Want to reorder terms in the input
      #   to match what's in the model, make sure the new documents have 0's for missing terms,
      #   then initialize the new terms.
      
      # 1. Copy the old term components
      beta <- list(beta=lapply(model$beta$logbeta, exp))
      if(!is.null(model$beta$kappa)) beta$kappa <- model$beta$kappa
      # 2. Reorder inputs to match order in model, putting new terms at the end
      # 3. Initialize new terms
    } else {
      cat("Vocabulary identical\n")
      beta <- list(beta=lapply(model$beta$logbeta, exp))
      if(!is.null(model$beta$kappa)) beta$kappa <- model$beta$kappa
    }
    # If there were new documents, initialize the entries in the document matrices
    #   now that we have updated term data.
    
    
    if(verbose) cat("Restarting Model...\n")
    #extract from a standard STM object so we can simply continue.

    
    # beta relates topics to terms

    sigma <- model$sigma # seems to be topic * topic matrix
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
  if (! is.null(mu$mu)) colnames(mu$mu) <- names(documents)
  if (! is.null(beta$kappa$m)) names(beta$kappa$m) <- vocab
  for (i in 1:length(beta$logbeta)) {
    colnames(beta$logbeta[[i]]) <- vocab
  }
  eta <- lambda[,-ncol(lambda), drop=FALSE]
  if (! is.null(eta)) rownames(eta) <- names(documents)
  theta <- exp(lambda - row.lse(lambda))
  if (! is.null(theta)) rownames(theta) <- names(documents)
  model <- list(mu=mu, sigma=sigma, beta=beta, settings=settings,
                vocab=vocab, convergence=convergence, 
                theta=theta, 
                eta=eta,
                invsigma=solve(sigma), time=time)
  class(model) <- "STM"  
  return(model)
}



    


