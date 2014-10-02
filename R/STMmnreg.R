#Penalized multinomial regression methods for content covariates
# using the distributed poisson approach

mnreg <- function(beta.ss,settings) {
  #Parse Arguments
  A <- settings$dim$A
  K <- settings$dim$K
  interact <- settings$kappa$interactions
  fixedintercept <- settings$kappa$fixedintercept
  alpha <- settings$tau$enet
  maxit <- settings$tau$maxit 
  nlambda <- settings$tau$nlambda
  lambda.min.ratio <- settings$tau$lambda.min.ratio
  ic.k <- settings$tau$ic.k
  thresh <- settings$tau$tol
  #Aggregate outcome data.
  counts <- do.call(rbind,beta.ss)

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
    covar <- sparseMatrix(veci, vecj, x=vecv)
  }
  
 
  if(fixedintercept) {  
    m <- settings$dim$wcounts$x
    m <- log(m) - log(sum(m))
  } else {
    m <- NULL #have to assign this to null to keep code simpler below
  }
  
  
  mult.nobs <- rowSums(counts) #number of multinomial draws in the sample
  offset <- log(mult.nobs)
  counts <- split(counts, col(counts))
  
  #########
  #Distributed Poissons
  #########
  
  #methods dispatch for S4 is crazy expensive so let's first define a function
  #for quickly extracting the coefficients from the model.
  subM <- function(x, p) {
    ind <- (x@p[p]+1):x@p[p+1]
    rn <- x@i[ind]+1
    y <- x@x[ind]
    out <- rep(0, length=nrow(x))
    out[rn] <- y
    out
  }
  
  #now do some setup of infrastructure
  verbose <- settings$verbose
  ctevery <- ifelse(length(counts)>100, floor(length(counts)/100), 1)
  out <- vector(mode="list", length=length(counts))
  #now iterate over the vocabulary
  for(i in 1:length(counts)) {
    #add fixed intercept if necessary
    if(is.null(m)) {
      offset2  <- offset
    } else {
      offset2 <- m[i] + offset    
    }
    #keep retrying glmnet until it works
    mod <- NULL
    while(is.null(mod)) {
      mod <- tryCatch(glmnet(x=covar, y=counts[[i]], family="poisson", 
                             offset=offset2, standardize=FALSE,
                             intercept=is.null(m), 
                             lambda.min.ratio=lambda.min.ratio,
                             nlambda=nlambda, alpha=alpha,
                             maxit=maxit, thresh=thresh),
        warning=function(w) return(NULL),
        error=function(e) stop(e))
      #if it didn't converge, increase nlambda paths by 20% 
      if(is.null(mod)) nlambda <- nlambda + floor(.2*nlambda)
    }
    dev <- (1-mod$dev.ratio)*mod$nulldev
    ic <- dev + ic.k*mod$df
    lambda <- which.min(ic)
    coef <- subM(mod$beta,lambda) #return coefficients
    if(is.null(m)) coef <- c(mod$a0[lambda], coef)
    out[[i]] <- coef
    if(verbose && i%%ctevery==0) cat(".")
  }
  if(verbose) cat("\n") #add a line break for the next message.  
  coef <- do.call(cbind, out)
  if(!fixedintercept) {
    #if we estimated the intercept add it in
    m <- coef[1,]
    coef <- coef[-1,]
  }
  kappa <- split(coef, row(coef))
  ##
  #predictions 
  ##
  #linear predictor
  linpred <- as.matrix(covar%*%coef)
  linpred <- sweep(linpred, 2, STATS=m, FUN="+")

  #softmax
  explinpred <- exp(linpred)
  beta <- explinpred/rowSums(explinpred)
  
  #wrangle into the list structure
  beta <- split(beta, rep(1:A, each=K))
  beta <- lapply(beta, matrix, nrow=K)
   
  kappa <- list(m=m, params=kappa)
  out <- list(beta=beta, kappa=kappa, nlambda=nlambda)
  return(out)
}
