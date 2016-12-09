selectModel <- function(documents, vocab, K,
                        prevalence, content, data=NULL,
                        max.em.its=100, verbose=TRUE, init.type = "LDA",
                        emtol= 1e-05, seed=NULL,runs=50, frexw=.7, 
                        net.max.em.its=2, netverbose=FALSE, M=10, N=NULL,
                        to.disk=F, ...){
  if(!is.null(seed)) set.seed(seed)

  if(is.null(N)){
    N <-  round(.2*runs)
  }
  
  if(runs<2){
    stop("Number of runs must be two or greater.")
  }
  
  if(runs<N){
    stop("Number in the net must be greater or equal to the number of final models.")
  }
  
  seedout <- NULL
  likelihood <- NULL
  cat("Casting net \n")
  for(i in 1:runs){
    cat(paste(i, "models in net \n"))
    mod.out <- stm(documents, vocab, K,
                   prevalence=prevalence, content=content, data=data, init.type=init.type,
                   max.em.its=net.max.em.its, emtol=emtol, verbose=netverbose,...)
    seedout[i] <- mod.out$settings$seed
    likelihood[i] <- mod.out$convergence$bound[length(mod.out$convergence$bound)]
  }

  keep <- order(likelihood, decreasing=T)[1:N]
  keepseed <- seedout[keep]
  cat("Running select models \n")
  runout <- list()
  semcoh <- list()
  exclusivity <- list()
  sparsity <- list()
  for(i in 1:length(keepseed)){
    cat(paste(i, "select model run \n"))
    initseed <- keepseed[i]
    mod.out <- stm(documents, vocab, K,
                   prevalence=prevalence, content=content, data=data, init.type=init.type, seed=initseed,
                   max.em.its=max.em.its, emtol=emtol, verbose=verbose,...)
    runout[[i]] <- mod.out
    if(to.disk==T){
      mod <- mod.out
      save(mod, file=paste("runout", i, ".RData", sep=""))
    }
    semcoh[[i]] <- semanticCoherence(mod.out, documents, M)
    if(length(mod.out$beta$logbeta)<2){
      exclusivity[[i]] <- exclusivity(mod.out, M=M, frexw=.7)
      sparsity[[i]] = "Sparsity not calculated for models without content covariates"
    }
    if(length(mod.out$beta$logbeta)>1){
      exclusivity[[i]] = "Exclusivity not calculated for models with content covariates"
      kappas <- t(matrix(unlist(mod.out$beta$kappa$params), ncol=length(mod.out$beta$kappa$params)))
      topics <-mod.out$settings$dim$K
      numsparse = apply(kappas[(K+1):nrow(kappas),], 1,function (x) sum(x<emtol))
      sparsity[[i]] = numsparse/ncol(kappas)
    }
  }
  return(list(runout=runout, semcoh=semcoh, exclusivity=exclusivity, sparsity=sparsity))
}
