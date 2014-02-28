#Penalized multinomial regression methods for content covariates
#(currently implements penalties between L1 and L2 will eventually
# incorporate between L0 and L1 with gamlr)

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
    m <- matrix(m, nrow=nrow(counts), ncol=ncol(counts), byrow=TRUE)
  } else {
    m <- NULL #have to assign this to null to keep code simpler below
  }
  
  mod <- glmnet(covar, counts, family="multinomial", intercept=!fixedintercept, 
                offset=m, maxit=maxit, standardize=FALSE, 
                alpha=alpha, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio,
                thresh=thresh)
  unpack <- unpack.glmnet(mod, nobs=nrow(counts), ic.k)
  coef <- unpack$coef
  kappa <- split(coef, row(coef))
  ##
  #predictions 
  ##
  #linear predictor
  linpred <- as.matrix(covar%*%coef)
  #if we estimated the intercept add it in  
  if(!fixedintercept) m <- unpack$intercept
  linpred <- m + linpred

  #softmax
  logbeta <- linpred - row.lse(linpred) 
  #wrangle into the list structure
  logbeta <- split(logbeta, rep(1:A, each=K))
  logbeta <- lapply(logbeta, matrix, nrow=K)
   
  kappa <- list(m=m[1,], params=kappa)
  out <- list(logbeta=logbeta, kappa=kappa)
  return(out)
}

#this is a small utility function to process glmnet results
#it will:
# 1) select the tuning parameter based on the information criterion weighted by k
#    (not in the sense of number of topics, in the sense of AIC k=2, BIC k=log n)
# 2) extract the coefficients.
# note that this would be straightforward with the coef() option but methods dispatch
# is way too slow with the S4 class.  This is particularly a problem for kappa with 
# large V
unpack.glmnet <- function(mod, nobs, ic.k) {
  dev <- (1-mod$dev.ratio)*mod$nulldev
  df <- colSums(mod$dfmat)
  ic <- dev + ic.k*df
  lambda <- which.min(ic)
  
  #methods dispatch here is crazy, so define out own function
  subM <- function(x, p) {
    ind <- (x@p[p]+1):x@p[p+1]
    rn <- x@i[ind]+1
    y <- x@x[ind]
    out <- rep(0, length=nrow(x))
    out[rn] <- y
    out
  }
  coef <- lapply(mod$beta, subM, lambda) #grab non-zero coefs
  coef <- do.call(cbind,coef)
  intercept <- mod$a0[,lambda]
  return(list(coef=coef, intercept=intercept))
}


