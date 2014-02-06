#######
# Functions for Inference of a Document in Logistic Normal Setting
#######
# File Contents
# 1) Doc Inference while temporarily collapsing over z
# 2) Likelihood/Gradient for Eta (Document-Topic Proportions)
# 3) Hessian for Joint Eta with Phi solved for within that set.
# 4) Approximate Bound

# Document Inference in Logistic-Normal Inference Collapsing Over z.
inferdoc.logisticnormal <- function(doc, priors, init) {
  doc.ct <- doc[2,]
  optim.out <- nlminb(start=init, objective=eta.likelihoodjoint, gradient=eta.gradientjoint,
                      doc.ct=doc.ct, mu=priors$mu,
                      siginv=priors$siginv, logbeta=priors$logbeta, Ndoc=sum(doc.ct))  
  
  #Solve for Hessian (to get nu)
  #   NB: this function actually computes a few other things to avoid the computations later including
  #   theta and phi.
  hessian <- eta.hessianjoint(optim.out$par, 
                              doc.ct=doc.ct, mu=priors$mu,
                              siginv=priors$siginv, logbeta=priors$logbeta, Ndoc=sum(doc.ct))
  
  #Now just putting things together in objects.  Note that we already negated the hessian
  eta <- list(lambda=optim.out$par, nu=solve(hessian$hessian)) 
  #Solve for the bound
  bound <- calc.bound(logbeta=priors$logbeta, eta=eta, theta=hessian$theta, priors=priors,doc.ct=doc.ct)
  #return all the document-level parameters.
  return(list(phis=hessian$phi, eta=eta, bound=bound)) 
}

##
# Functions for Joint Eta Inference
##
# Note the functions below are quite fast but are slightly less readable than
# archived versions.  Extensive testing was done for the gradient.  Gains
# could be made by caching calculations between the two but that's quite a bit harder.
# We implemented earlier slightly less efficient versions of this in C++ and didn't observe
# any speed gains, although it may be worth retesting at some point with new functions.

eta.likelihoodjoint <- function(eta, doc.ct, mu, siginv,logbeta, Ndoc=sum(doc.ct)) {
  eta.0 <- c(eta,0)
  # {\sum_{v=1}^V c_v log [\sum_k beta_{k,v} exp(eta_k)] }- Wlog \sum exp(eta_k)
  lse.piece <- col.lse(eta.0 +logbeta)
  part1.1 <- sum(doc.ct*lse.piece)
  part1.2 <- Ndoc*lse(eta.0)
  # -1/2 (eta - mu)^T Sigma (eta - mu)
  part2 <- -.5*(eta-mu)%*%siginv%*%(eta-mu) 
  #note that above we remove transpose from initial eta-mu in order to get some minor speed boosts.
  -(part1.1 - part1.2 + part2)  
} 

#Gradient of the Joint Objective
# This function is called a lot so I try to heavily optimize it here.  
# note that there are several ways to do the normalization in the fourth line
# I experimented with versions with and without transpose, sweep, and rep(denom,each=)
# style versions..  The transpose version is the fastest presumably because the matrices
# are small and colSums() is marginally faster than rowSums()
# Presumably these calculations could be cached between the objective and gradient but that's
# a more complicated set of work.
eta.gradientjoint <- function(eta, doc.ct, mu, siginv,logbeta, Ndoc=sum(doc.ct)) {
  eta.0 <- c(eta, 0)
  Ez <- exp(eta.0 + logbeta) 
  denom <- colSums(Ez)/doc.ct
  Ez <- t(Ez)/denom
  part1 <- colSums(Ez) - Ndoc*softmax(eta.0)
  part2 <- - siginv%*%(eta-mu) 
  -(part1[-length(eta.0)] + part2)
}

#This function calculates the Hessian and phi
eta.hessianjoint <- function(eta,doc.ct, mu, siginv, logbeta, Ndoc=sum(doc.ct)) {
  #basic transforms
  eta.0 <- c(eta, 0)
  theta <- softmax(c(eta.0))
  
  #pieces for the derivatives of the exp(eta)beta part
  EB <- exp(eta.0 + logbeta) #calculate exp(eta)\beta for each word
  EB <- t(EB)/colSums(EB) #transpose and norm by (now) the row
  #at this point EB is the phi matrix
  phi <- EB*(doc.ct) #multiply through by word count
  phisums <- colSums(phi)
  phi <- t(phi) #transpose so its in the K by W format expected
  EB <- EB*sqrt(doc.ct) #set up matrix to take the cross product
  
  #First piece is the quotient rule portion that shows up from E[z], second piece is the part
  # that shows up regardless as in Wang and Blei (2013) for example.  Last element is just siginv
  hess <- -((diag(phisums) - crossprod(EB)) -Ndoc*(diag(theta) - theta%o%theta))[1:length(eta),1:length(eta)] + siginv
  return(list(hessian=hess, phi=phi, theta=theta))
}

calc.bound <- function(logbeta, eta, theta, priors,doc.ct) {
  logphinorm <- col.lse(log(theta) + logbeta)
  part1 <- sum(doc.ct*logphinorm)

  likelihood <- part1 + 
    .5*determinant(eta$nu, logarithm=TRUE)$modulus -
    .5*t(eta$lambda-priors$mu)%*%priors$siginv%*%(eta$lambda-priors$mu) -
    priors$sigmaentropy
  as.numeric(likelihood)
}

