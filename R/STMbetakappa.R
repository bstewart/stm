
#This function controls all the topic estimation functions.
opt.beta <- function(beta.ss, kappa=NULL, settings=NULL) {
  if(is.null(kappa)) {
    return(list(logbeta=list(opt.betaLDA(beta.ss[[1]]))))
  } else {
    if(settings$tau$mode=="L1") {
      return(mnreg(beta.ss,settings))
    } else {
      return(estKappa(beta.ss, kappa, settings))
    }
  }
} 

#Estimation of beta in LDA style.
opt.betaLDA <- function(beta.ss) {
  rowTotals <- rowSums(beta.ss)
  safelog(beta.ss) - log(rowTotals)
}