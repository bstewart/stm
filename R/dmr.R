# Dirichlet Multinomial Regression model
# (Workhorse function, also pretty fast.)

#NB: I've profiled this a few times now.  Its about as fast as its going to be in R
# because the cost is in the indexing of components.  I tested a few different data structures
# but it didn't really make a difference.

dmr.control <- function(documents,vocab, K, X, 
                        sig2=c(10, rep(.5, times=(ncol(X)-1))), 
                        eta=.01, 
                        tol=.01,maxits=10000,
                        verbose=TRUE, printevery=5) {
  
  #Create the data.
  # three equal length vectors: the document id, word indices, the count of each word
  # note this is the simple-triplet-matrix i,j,v respectively.
  dat <- doc.to.ijv(documents)
  
  #test for -indexing
  if(0 %in% dat$j) {
    ldadoc <- ijv.to.doc(dat$i,dat$j,dat$v)
  } else {
    ldadoc <- ijv.to.doc(dat$i,dat$j-1,dat$v)
  }
  
  ##
  #Initializing
  ##
  #Need to track four things:
  # phi  is the K by V matrix of counts of words by topic
  # theta is the D length list of expected counts in each document
  # NZ is the is the K length vector containing counts in the topics.
  # gamma is the token length list, where each element is a K-length vector
  
  #Use a few passes of Gibbs to get things started
  init <- lda.collapsed.gibbs.sampler(ldadoc,K,vocab, num.iterations=100, alpha=.01, eta=.01)
  #NB: no point in doing burnin because the statistics have to sink up and it doesn't give us the 
  #    proper values for phi, so it isn't possible to work back a meaningful use of those.
  
  #theta
  x <- init$document_sums
  theta <- split(x, col(x))
  rm(x)
  #phi
  phi <- init$topics
  
  #topic sums
  NZ <- rowSums(phi)
  #gamma
  gvec <- unlist(init$assignments)+1 #offset the zero index.
  gvec <- split(gvec, 1:length(gvec))
  gamma <- lapply(gvec, function(x) {
                  out <- rep(0,K)
                  out[x] <- 1
                  out
                  })
  rm(gvec)
  
  #initialize the priors
  alpha <- rep(list(rep(.01,K)), length(documents)) #list of document alphas
  etasum <- length(phi)*eta
  
  ###
  # Define Optimization Functions for Later Use
  ###
  
  lambda.likelihood <- function(lambda, X, sig2, thetamat, Nd) {
    #Note it seems to be pretty sensitive to initialization.
    #take lambda and map from a vector to block form
    lambda <- matrix(lambda, ncol=ncol(thetamat))
    fit <- exp(X%*%lambda)
    dSums <- rowSums(fit)
    part1 <- sum(lgamma(dSums)) - sum(lgamma(dSums + Nd))
    part2 <- sum(lgamma(fit + thetamat) - lgamma(fit)) 
    part3 <-  sum(-lambda^2/(2*sig2)) - ncol(thetamat)*sum(log(sqrt(2*pi*sig2)))
    part1 + part2 + part3
  }
  
  lambda.gradient <- function(lambda, X, sig2, thetamat, Nd) {
    #take lambda and map from a vector to block form
    #Line numbers below reference equation 2 in Mimno and McCallum
    lambda <- matrix(lambda, ncol=ncol(thetamat))
    fit <- exp(X%*%lambda)
    dSums <- rowSums(fit)
    part2 <- digamma(dSums) - digamma(dSums + Nd) #line 2 in DMR Paper
    
    part3 <- digamma(fit + thetamat) - digamma(fit) #line 2 + line 3
    part1 <- fit*(part2 + part3) #all by the first term of line 1, multiplied by lines 2 and 3
    out <- t(X)%*%part1  #all pieces except the prior
    out <- out - lambda/sig2 #including the prior on the coefficients
    return(as.numeric(out))
  }
  
  ###
  # Perform Inference
  ###
  # Intuition:
  # the gamma update is two ratios multiplied:
  # first is probability of w_i under topic j (thus need word/topic count matrix and topic-count matrix)
  # second is probability of topic j in document d. (thus needing expected topic count in doc)
  
  # Steps:
  # At each word index:
  #  1) pull the column of phi
  #  2) pull the list element of theta
  #  3) pull the previous gamma value.
  #  4) Adjust: phi, theta and NZ by count and previous value of gamma
  #  5) Calculate gamma
  #  6) update ss
  #  7) update gamma
  
  # Notes:
  # we are just using a for loop to avoid having to do global reassignment out of a function.
  
  #We start by defining a bunch of useful objects/constants
  convergence <- c(100)
  n.passes <- 1
  nwords <- length(gamma) #number of unique document-type combinations
  docids <- dat$i
  windex <- dat$j
  counts <- dat$v
  docsums <- dat$rowsums
  rm(dat)
  
  #Now begin loop
  while(convergence[length(convergence)] > tol & n.passes <=maxits) {
    #save the phis for a convergence comparison later.
    compare.phi <- phi
    
    ##
    # E-Step
    for (i in 1:nwords) {
      #grab statistics corresponding to the word
      #(using convention of including w in front of full name to indicate "working" version)
      wordindex <- windex[i]
      docid <- docids[i]
      count <- counts[i]
      
      wphi <- phi[,wordindex]
      wtheta <- theta[[docid]]
      wgamma <- gamma[[i]]
      docsum <- docsums[docid] #number of words in the document
      docalpha <- alpha[[docid]]
        
      #adjust statistics to remove influence of current word
      remove <- wgamma*count
      wphi <- wphi - remove
      wtheta <- wtheta - remove
      NZ <- NZ - remove
      docsum <- docsum - count
      
      #update gamma
      p1 <- (wphi + eta)/(NZ + etasum)
      p2 <- (wtheta + docalpha)/(docsum + docalpha*K)        
      wgamma <- p1*p2
      wgamma <- wgamma/sum(wgamma)
      
      #update global records
      update <- wgamma*count
      phi[,wordindex] <- wphi + update
      theta[[docid]] <- wtheta + update
      NZ <- NZ + update
      gamma[[i]] <- wgamma
    }
    #convergence check
    newconv <- sum(abs(phi/rowSums(phi) - compare.phi/rowSums(compare.phi))/K)
    convergence <- c(convergence,newconv)
 
    ##
    # M-Step
    thetamat <- do.call(rbind,theta)
    init <- rep(0, ncol(X)*K)
    optim.out <- optim(par=init, fn=lambda.likelihood,
                       gr=lambda.gradient, method="L-BFGS-B",
                       control=list(fnscale=-1, maxit=1000),
                       X=X, sig2=sig2, thetamat=thetamat, Nd=docsums)
   
    lambda <- matrix(optim.out$par, ncol=ncol(thetamat))
    alpha <- exp(X%*%lambda)
    alpha <- split(alpha,row(alpha))
    rm(thetamat)
   
    if(verbose && (n.passes%%printevery==0)) {
      cat(paste("Finishing Iteration ",n.passes, "\n"))
    }
    n.passes <- n.passes + 1
  }
  
  ###
  # Package and return
  ###
  cat(paste("Finished After ", n.passes, " iterations."))
  phi <- phi/rowSums(phi)
  theta <- do.call(rbind,theta)
  theta <- theta/rowSums(theta)
  return(list(gamma=gamma, topic_sums=NZ, topics=phi, document_sums=theta, 
              convergence=convergence, coefs=lambda, alpha=alpha, 
              vocab=vocab))
}



