#Initializing Functions
stm.initialize <- function(documents, settings) {
  beta.init <- NULL 
  K <- settings$dim$K
  V <- settings$dim$V
  A <- settings$dim$A
  N <- settings$dim$N
  mode <- settings$init$mode
  #Special Beta Initializations
  if(mode=="LDA") {
      beta.init <- init.params.LDA(documents,K, V, 
                                   alpha=1,
                                   eta=.1)
  }
  if(mode=="DMR") {
    dmr.mod <-  dmr.control(documents,settings$vocab, K, settings$covariates$X,maxits=100,verbose=FALSE)
    beta.init <- log(dmr.mod$phi)
  }
  if(mode=="User") {
    beta.init <- settings$init$userinit
    if(is.matrix(beta.init)) beta.init <- list(beta.init)
  }

  
  model <- init.params(K, N, V, A, beta.init)    
 
  if(!settings$kappa$LDAbeta) {
    model$beta$kappa <- kappa.init(documents, K, V, A, interactions=settings$kappa$interactions)
  }
  model$settings <- settings
  return(model)
}

#Initialize the parameters 
init.params <- function(K, N, V, A, logbeta=NULL) {
  
  mu <- matrix(0, nrow=(K-1),ncol=1)
  sigma <- diag(20, nrow=(K-1))

  if(is.null(logbeta)) {
    logbeta <- vector(mode="list",length=A)
    for (i in 1:A) {
      beta <- matrix(rgamma(V * K, .1), ncol = V)
      logbeta[[i]] <- log(beta) - log(rowSums(beta))
    }
  } else{
    if(length(logbeta)!=A) {
      if(!is.list(logbeta)) logbeta <- list(logbeta)
      logbeta <- rep(logbeta,A)    
    }
  }
  out <- list(mu=list(mu=mu), sigma=sigma, beta=list(logbeta=logbeta))
  return(out)
}

#Initialize the parameters using a few passes of collapsed Gibbs Sampling. 
init.params.LDA <- function(documents, K, V, alpha, eta) {
  resetdocs <- documents
  for (i in 1:length(documents)) {
    resetdocs[[i]][1,] <- as.integer(resetdocs[[i]][1,]-1)
  }
  
  mod <- lda.collapsed.gibbs.sampler(resetdocs, K, 1:V, num.iterations=50, 
                                     alpha=alpha, eta=eta)
  
  #smoothing a bit to avoid sharp distributions
  topics <- as.numeric(mod$topics) + 1 #note the recasting here is to avoid integer overflow on really large doc sets
  topics <- matrix(topics, nrow=K)
  logbeta <- log(topics) - log(rowSums(topics))
  return(list(logbeta))
}

###
# Kappa initialization
###
kappa.init <- function(documents, K, V, A, interactions) {
  kappa.out <- list()
  #Calculate the baseline log-probability (m)
  freq <- matrix(unlist(documents),nrow=2) #break it into a matrix
  freq <- split(freq[2,], freq[1,]) #shift into list by word type
  m <- unlist(lapply(freq, sum)) #sum over the word types
  m <- m/sum(m)
  #m <- log(m)
  m <- log(m) - log(mean(m)) #logit of m
  kappa.out$m <- m
  
  #Defining parameters
  aspectmod <- A > 1
  if(aspectmod) {
    interact <- interactions 
  } else {
    interact <- FALSE
  }
  
  #Create the parameters object
  parLength <- K + A*aspectmod + (K*A)*interact
  kappa.out$params <- vector(mode="list",length=parLength)
  for(i in 1:length(kappa.out$params)) {
    kappa.out$params[[i]] <- rep(0, V)
  }
  
  #Create a running sum of the kappa parameters starting with m
  kappa.out$kappasum <- vector(mode="list", length=A)
  for (a in 1:A) {
    kappa.out$kappasum[[a]] <- matrix(m, nrow=K, ncol=V, byrow=TRUE)
  }
  
  #create covariates. one element per item in parameter list.
    #generation by type because its conceptually simpler
  if(!aspectmod & !interact) {
    kappa.out$covar <- list(k=1:K, a=rep(NA, parLength), type=rep(1,K))
  }
  if(aspectmod & !interact) {
    kappa.out$covar <- list(k=c(1:K,rep(NA,A)), a=c(rep(NA, K), 1:A), type=c(rep(1,K), rep(2,A)))      
  }
  if(interact) {
    kappa.out$covar <- list(k=c(1:K,rep(NA,A), rep(1:K,A)), 
                        a=c(rep(NA, K), 1:A, rep(1:A,each=K)), 
                        type=c(rep(1,K), rep(2,A), rep(3,K*A)))            
  }
  return(kappa.out)
}