####
# Fast document level inference
####
#these functions are considerably faster than their predecessors as a
#consequence of (a) integrating out z, (b) avoiding expensive log-sum-exp calls.

# Document Inference with Logistic-Normal
#  (assumes that beta has already been subset)
logisticnormal <- function(eta, mu, siginv, beta, doc, sigmaentropy) {
  doc.ct <- doc[2,]
  Ndoc <- sum(doc.ct)
  #even at K=100, BFGS is faster than L-BFGS
  optim.out <- optim(par=eta, fn=lhood, gr=grad,
                     method="BFGS", control=list(maxit=500),
                     doc.ct=doc.ct, mu=mu,
                     siginv=siginv, beta=beta, Ndoc=Ndoc)
  
  #Solve for Hessian/Phi/Bound returning the result
  hpb(optim.out$par, doc.ct=doc.ct, mu=mu,
      siginv=siginv, beta=beta, Ndoc=Ndoc,
      sigmaentropy=sigmaentropy)
}

##
# Fast likelihood calculation
##
#fairly carefully engineered
#benchmarking notes: 
# 1) Part 1
#    Some brief profiling done to test that the formulation colSums(beta*expeta) was better
#    than the matrix multiplication. (NB: I returned to this a second time without reading
#    this note.  for K=100 this is still true although difference is small.)
#    Benchmarking suggests that the expense here is in the multiplication, with the only other
#    non-negligible part being the colSums.  This suggests caching the multiplication to be shared
#    with gradient could be useful.
# 2) Part 2
#    for the part 2 piece the option (-.5*diff%*%siginv%*%diff) is significantly faster 
#    for small samples but at K=50 its about 3 times faster to use the 
#    formulation with the crossproduct.  I in turn found that for K=100 it was better
#    to use the matrix multiply.  I tested a verison with a Cholesky decomp but it didn't help.
lhood <- function(eta, doc.ct, mu, siginv,beta, Ndoc=sum(doc.ct)) {
  # {\sum_{v=1}^V c_v log [\sum_k beta_{k,v} exp(eta_k)] }- Wlog \sum exp(eta_k)
  expeta <- c(exp(eta),1)
  part1 <- sum(doc.ct*log(.colSums(beta*expeta, nrow(beta), ncol(beta)))) - Ndoc*log(sum(expeta))
  # -1/2 (eta - mu)^T Sigma (eta - mu)
  diff <- eta-mu
  part2 <- .5*sum(diff*(siginv %*% diff))
  part2 - part1  
} 

#Gradient of the Joint Objective
#  again fairly carefully engineered.  
#  - Ez reused to minimize memory
#  - the denominator is flipped so we can multiply which is faster.
#  - column is dropped off the first time it is computationally advantageous
#  - second component tested with colSums and * but its much faster as a matrix multiply
#  - we replaced a version with the transpose for the more efficient multiplication.
grad <- function(eta, doc.ct, mu, siginv, beta, Ndoc=sum(doc.ct)) {
  expeta.sh <- exp(eta) 
  expeta <- c(expeta.sh,1)
  Ez <- expeta*beta
  denom <- doc.ct/.colSums(Ez, nrow(Ez), ncol(Ez))
  part1 <- (Ez%*%denom)[-length(expeta)] - expeta.sh*(Ndoc/sum(expeta))  
  part2 <- siginv%*%(eta-mu) 
  as.numeric(part2 - part1)
}

# Hessian/Phi/Bound
#   NB: Hessian function is not as carefully benchmarked as it isn't called
#       nearly as often.  Particularly I suspect some of the elements in the
#       cross product could be sped up quite a bit.
#   NB: Bound and hessian barely communicate with one another here.
#       1/29/15 I revisited this and tuned it up a bit but was unable to 
#       substantially reduce the communication issues.
#   NB: 1/29/15- benchmarking on poliblog5k with 100 topics suggests that
#       crossprod() is a surpringly big chunk of the cost despite it only
#       appearing in the Hessian.  Its the single bigest chunk after "*"
#       its producing a K by K matrix but still surprising...
hpb <- function(eta, doc.ct, mu, siginv, beta, Ndoc=sum(doc.ct), sigmaentropy) {
  #basic transforms
  expeta <- c(exp(eta),1)
  theta <- expeta/sum(expeta)
  
  #pieces for the derivatives of the exp(eta)beta part
  EB <- expeta*beta #calculate exp(eta)\beta for each word
  EB <- t(EB)/colSums(EB) #transpose and norm by (now) the row
  
  #at this point EB is the phi matrix
  phi <- EB*(doc.ct) #multiply through by word count
  phisums <- colSums(phi)
  phi <- t(phi) #transpose so its in the K by W format expected
  EB <- EB*sqrt(doc.ct) #set up matrix to take the cross product
  
  Ntheta <- Ndoc*theta
  rtNtheta <- sqrt(Ndoc)*theta
  #First piece is the quotient rule portion that shows up from E[z], second piece is the part
  # that shows up regardless as in Wang and Blei (2013) for example.  Last element is just siginv
  hess <- -((diag(phisums) - crossprod(EB)) - 
              (diag(Ntheta) - rtNtheta%o%rtNtheta))[1:length(eta),1:length(eta)] + siginv
  
  ###
  # Bound
  
  cobj <- try(chol.default(hess), silent=TRUE)
  if(class(cobj)=="try-error") {
    #brute force solve
    nu <- solve(hess)
    #only if we would produce negative variances do we bother doing nearPD
    if(any(diag(nu)<0)) nu <- as.matrix(nearPD(nu)$mat)
    detTerm <- .5*determinant(nu, logarithm=TRUE)$modulus
  } else {
    #this is the normal case where the hessian an be decomposed
    nu <- chol2inv(cobj)
    #this is equivalent to the 1/2 of the log determinant as given above
    detTerm <- -sum(log(diag(cobj)))
  }
  diff <- eta - mu
  
  logphinorm <- log(as.numeric(theta%*%beta))
  part1 <- sum(doc.ct*logphinorm)
  bound <- part1 + detTerm - #the determinant is now pre-computed so we can speed it up where possible.
           .5*sum(diff*crossprod(diff,siginv)) -
           sigmaentropy
  bound <- as.numeric(bound)
    
  #bundle everything up.
  return(list(phis=phi, eta=list(lambda=eta, nu=nu), bound=bound))
}
