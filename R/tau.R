
opt.tau <- function(kappa.init=NULL, i=NULL, kappa=NULL, mode, fixedprior=NULL) {
  if(mode=="Jeffreys") {
    return(jeffreys(kappa.init))
  }
  if(mode=="Pooled") {
    return(pooledvar(kappa,i))
  }
  if(mode=="Fixed") {
    paramgroup <- which(kappa$covar$type==kappa$covar$type[i])
    return(fixedprior[paramgroup])
  }
}

##
#This method implements a simple pooled MAP estimate by type
##
pooledvar <- function(kappa, i) {
  paramgroup <- which(kappa$covar$type==kappa$covar$type[i])
  tau <- mean(as.numeric(unlist(kappa$params[paramgroup]))^2)
  prec <- 1/tau
  prec <- ifelse(prec > 1e5,1e5,prec)
  return(prec)
}

###
# Function to estimate jeffreys variance.
###
jeffreys <- function(x) {
  if(sum(abs(x))==0) {
    return(1)
  } else {
    x <- x^2
    prec <- 1/x
    prec[prec>1e5] <- 1e5
    return(prec)
  }
}