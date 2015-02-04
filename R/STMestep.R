#E-Step for a Document Block
#[a relatively straightforward rewrite of previous
# code with a focus on avoiding unnecessary computation.]

#Input: Documents and Key Global Parameters
#Output: Sufficient Statistics

# Approach:
# First we pre-allocate memory, and precalculate where possible.
# Then for each document we:
#  (1) get document-specific priors, 
#  (2) infer doc parameters, 
#  (3) update global sufficient statistics
# Then the sufficient statistics are returned.

#Let's start by assuming its one beta and we may have arbitrarily subset the number of docs.
estep <- function(documents, beta.index, update.mu, #null allows for intercept only model  
                       beta, lambda.old, mu, sigma, 
                       verbose) {
  
  #quickly define useful constants
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  ctevery <- ifelse(N>100, floor(N/100), 1)
  if(!update.mu) mu.i <- as.numeric(mu)
  
  # 1) Initialize Sufficient Statistics 
  sigma.ss <- diag(0, nrow=(K-1))
  beta.ss <- vector(mode="list", length=A)
  for(i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
  
  # 2) Precalculate common components
  sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
  siginv <- solve(sigma)
    
  # 3) Document Scheduling
  # For right now we are just doing everything in serial.
  # the challenge with multicore is efficient scheduling while
  # maintaining a small dimension for the sufficient statistics.
  print("Using the new code")
  for(i in 1:N) {
    #update components
    doc <- documents[[i]]
    words <- doc[1,]
    aspect <- beta.index[i]
    init <- lambda.old[i,]
    if(update.mu) mu.i <- mu[,i]
    beta.i <- beta[[aspect]][,words,drop=FALSE]
    
    #infer the document - formerly the logisticnormal function
    #even at K=100, BFGS is faster than L-BFGS
    expeta <- 0
    diff <- 0
    doc.ct <- doc[2,]
    Ndoc <- sum(doc.ct)
    cols = ncol(beta.i)
    rows = nrow(beta.i)
    siginvdiff <- 0
    betaexpeta_colsums <- 0
    betaexpeta <- 0
    sumexpeta <- 0
    global_eta <- 0
    eta <- optim(par=init, 
                       fn=function(eta) {
                         # {\sum_{v=1}^V c_v log [\sum_k beta_{k,v} exp(eta_k)] }- Wlog \sum exp(eta_k)
                         global_eta <<- eta
                         expeta <<- c(exp(eta),1)
                         sumexpeta <<- sum(expeta)
                         betaexpeta <<- beta.i * expeta
                         betaexpeta_colsums <<- .colSums(betaexpeta, 
                                                         rows, 
                                                         cols)
                         part1 <- sum(doc.ct*log(betaexpeta_colsums)) - Ndoc*log(sumexpeta)
                         # -1/2 (eta - mu)^T Sigma (eta - mu)
                         diff <<- eta-mu.i
                         siginvdiff <<- siginv %*% diff
                         part2 <- .5*sum(diff*siginvdiff)
                         part2 - part1  
                       },  gr=function(eta) {
                         if (identical(eta,global_eta)) {
                           denom <- doc.ct/betaexpeta_colsums
                           part1 <- (betaexpeta%*%denom)[-length(expeta)] - expeta[1]*(Ndoc/sumexpeta)  
                           as.numeric(siginvdiff - part1)
                         } else {
                           expeta.sh <- exp(eta) 
                           expeta <- c(expeta.sh,1)
                           Ez <- expeta*beta
                           denom <- doc.ct/.colSums(Ez, nrow(Ez), ncol(Ez))
                           part1 <- (Ez%*%denom)[-length(expeta)] - expeta.sh*(Ndoc/sum(expeta))  
                           part2 <- siginv%*%(eta-mu) 
                           as.numeric(part2 - part1)                           
                         }
                       },
                       method="BFGS", 
                       control=list(maxit=500)
    )$par
    
#     doc.results <- logisticnormal(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, 
#                                  doc=doc, sigmaentropy=sigmaentropy)
    
    #Solve for Hessian/Phi/Bound returning the result
      theta <- expeta/sumexpeta
  
   #pieces for the derivatives of the exp(eta)beta part
       EB <- t(betaexpeta)/betaexpeta_colsums #transpose and norm by (now) the row
  
   #at this point EB is the phi matrix
      phi <- EB*(doc.ct) #multiply through by word count
      phisums <- colSums(phi)
      phi <- t(phi) #transpose so its in the K by W format expected
      EB <- EB*sqrt(doc.ct) #set up matrix to take the cross product
    
    #First piece is the quotient rule portion that shows up from E[z], second piece is the part
    # that shows up regardless as in Wang and Blei (2013) for example.  Last element is just siginv
     hess <- -((diag(phisums) - crossprod(EB)) - 
                 Ndoc*(diag(theta) - theta%o%theta))[1:length(eta),1:length(eta)] + siginv
    
    ###
    # Bound
    
    nu <- try(chol2inv(chol.default(hess)), silent=TRUE)
    if(class(nu)=="try-error") {
      #brute force solve
      nu <- solve(hess)
      #only if we would produce negative variances do we bother doing nearPD
      if(any(diag(nu)<0)) nu <- as.matrix(nearPD(nu)$mat)
    }
     logphinorm <- log(colSums(theta*beta.i))
     part1 <- sum(doc.ct*logphinorm)

    bound[i] <- as.numeric(part1 + .5*determinant(nu, logarithm=TRUE)$modulus -
      .5*sum(diff*crossprod(diff,siginv)) - sigmaentropy)
    
    # done inferring the document

    # update sufficient statistics 
    sigma.ss <- sigma.ss + nu
    beta.ss[[aspect]][,words] <- phi + beta.ss[[aspect]][,words]
#     bound[i] <- doc.results$bound
    lambda[[i]] <- eta

    if(verbose && i%%ctevery==0) cat(".")
  }
  if(verbose) cat("\n") #add a line break for the next message.
  
  #4) Combine and Return Sufficient Statistics
  lambda <- do.call(rbind, lambda)
  return(list(sigma=sigma.ss, beta=beta.ss, bound=bound, lambda=lambda))
}
