#Gamma Lasso using the TextIR package
gammaLasso <- function(beta.ss,settings) {
  A <- settings$dim$A
  K <- settings$dim$K
  interact <- settings$kappa$interactions
  
  counts <- do.call(rbind,beta.ss)
  
  #Let's do a simple triplet matrix for the covariates
  #Three Cases
  if(A==1) { #Topic Model
    covar <- diag(1, nrow=K)
  }
  if(A!=1) { #Topic-Aspect Models
    #Topics
    veci <- 1:nrow(counts)
    vecj <- rep(1:K,A)
    #aspects
    veci <- c(veci,1:nrow(counts))
    vecj <- c(vecj,rep((K+1):(K+A), each=K))
    if(interact) {
      veci <- c(veci, 1:nrow(counts))
      vecj <- c(vecj, (K+A+1):(K+A+nrow(counts)))
    }
    vecv <- rep(1,length(veci))
    covar <- simple_triplet_matrix(veci, vecj, vecv)
  }
  
  #Okay now run the model
  mod <- suppressWarnings(dmr(cl=NULL, counts=counts, 
                              covars=covar,
                              normalize=FALSE, penalty=c(shape=3,rate=1))) 
  
  #Prep Results for Output
  logprob <- log(predict(mod,covar))
  logbeta <- vector(mode="list", length=A)
  indx <- rep(1:A, each=K)
  for(i in 1:A) {
    logbeta[[i]] <- logprob[indx==i,]
  }
  rm(logprob)
  
  kappa <- list()
  coef <- coef(mod)
  kappa$m <- coef[1,]
  kappa$params <- vector(mode="list", length=(nrow(coef)-1))
  for(i in 1:length(kappa$params)) {
    kappa$params[[i]] <- as.numeric(coef[(1+i),])
  }  
  out <- list(logbeta=logbeta, kappa=kappa)
  return(out)
}