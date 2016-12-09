#########################################################################
#         Generalized Model Stability Analysis Routines for STM         #
#             Antonio Coppola, Harvard University, July 2014            #
#########################################################################

multiSTM <- function(mod.out=NULL, ref.model=NULL, 
                     align.global=FALSE, mass.threshold=1, 
                     reg.formula=NULL, metadata=NULL, reg.nsims=100,
                     reg.parameter.index=2, verbose=TRUE,
                     from.disk=FALSE){
  # Catch bad inputs
  if(from.disk==F && is.null(mod.out))
    stop("Must supply a selectModel output object if not loading models from disk.")
  if(from.disk==T && !is.null(mod.out))
    stop("Loading models from disk and RAM are mutually exclusive options.")
  if(align.global && !requireNamespace("clue", quietly=TRUE)) 
    stop("Please install the clue package to use align.global=TRUE")
  if (findInterval(mass.threshold, c(0,1), rightmost.closed=TRUE) != 1)
    stop("Mass threshold must be in the [0,1] range.")
  if (!is.null(reg.formula) & is.null(metadata))
    stop("Must provide metadata in order to perform covariate effect analysis.")
  if(length(mod.out$runout[[1]]$beta$logbeta)>1) 
    stop("multiSTM does not yet work with content covariates.")
  
  # Load from disk if necessary
  if(from.disk==T){
    fnames <- list.files()
    datasets <- fnames[which(grepl("*.RData", fnames))]
    datasets <- datasets[order(datasets)]
    mod.out <- list()
    mod.out$runout <- list()
    for(i in seq(length(datasets))){
      mod <- load(datasets[i])
      mod.out$runout[[i]] <- get(mod)
    }
  }
  
  # Number of models in object
  N <- length(mod.out$runout)
  if(!is.null(ref.model)){
    if(ref.model < 1 || ref.model > N) 
      stop("Reference model number outside bounds.")
  } 
  
  # Grab the bounds
  lb <- unlist(lapply(mod.out$runout, function(x) max(x$convergence$bound)))
  
  # If no reference model is specified grab the max
  if(is.null(ref.model)) {
    glob.max <- which.max(lb)
  } else {
    glob.max <- ref.model #don't love this nomenclature but keeping it for consistency -bs
  }
  # Load reference objects
  ref.mod <- mod.out$runout[[glob.max]]
  ref.beta <- exp(ref.mod$beta$logbeta[[1]])
  ref.theta <- exp(ref.mod$theta)
  ref.thoughts <- findThoughts(ref.mod, n=10)
  ref.words <- apply(ref.beta, 1, order, decreasing=TRUE)[1:10,]
  
  # Load topic number
  K <- ref.mod$settings$dim$K
  
  # Setup for covariate effect stability analysis
  if(!is.null(reg.formula)){
    
    # Perfom reference regression
    ref.reg <- estimateEffect(reg.formula, ref.mod,
                              metadata=metadata,
                              uncertainty="Global",
                              nsims=reg.nsims)
    
    # Instantiate covariance matrix
    var.matrix <- matrix(NA, N, K)
    cov.effects <- list()
    length(cov.effects) <- N
    
    # Confidence intervals for regression estimates
    confidence.intervals <- matrix(NA, 2, K)
    for (j in 1:K){
      mean <- ref.reg$parameters[[j]][[1]]$est[reg.parameter.index]
      stdev <- sqrt(ref.reg$parameters[[j]][[1]]$vcov[reg.parameter.index,
                                                      reg.parameter.index])
      confidence.intervals[1, j] <- mean - 1.96*stdev
      confidence.intervals[2, j] <- mean + 1.96*stdev
    }
  }
  
  # Initialize objects of interest
  L1 <- list()
  thoughtct <- list()
  wordct <- list()
  
  # Main loop over models
  for (i in seq(N)){
    
    # Load beta matrix 
    betamod <-  exp(mod.out$runout[[i]]$beta$logbeta[[1]])
    if(ncol(betamod)!=ncol(ref.beta)) 
      next
    
    # Perform topic alignment
    if(align.global==TRUE)
      align <- clue::solve_LSAP(ref.beta%*%t(betamod), maximum=TRUE)
    else{
      product <- ref.beta%*%t(betamod)
      align <- apply(product,1,which.max)
    }
    
    # Compute L1
    L1[[i]] <- rowSums(abs(betamod[align,] - ref.beta))
    
    # Distance in words
    words <- apply(betamod, 1, order, decreasing=TRUE)[1:10,align]
    topwords <- vector(length=ncol(words))
    for(j in 1:ncol(words))
      topwords[j] <- sum(ref.words[,j] %in% words[,j])
    wordct[[i]] <- topwords
    
    # Shared top documents
    thoughts <- findThoughts(mod.out$runout[[i]], n=10)    
    thoughts <- thoughts[[1]][align]
    topdocs <- vector(length=length(thoughts))
    for(k in 1:length(thoughts))
      topdocs[k] <- sum(ref.thoughts[[1]][[k]] %in% thoughts[[k]])
    thoughtct[[i]] <- topdocs
    
    # Covariate effect stability analysis
    if (!is.null(reg.formula)){
      cov.effects[[i]] <- ref.reg <- estimateEffect(reg.formula, 
                                                    mod.out$runout[[i]],
                                                    metadata=metadata, 
                                                    uncertainty="Global",
                                                    nsims=reg.nsims)
      for(h in 1:K){
        var.matrix[i, h] <- cov.effects[[i]]$parameters[[h]][[1]]$est[reg.parameter.index]
      }
      var.matrix[i, ] <- var.matrix[i, align]
    }
    
    # Report progress
    if (verbose==TRUE & i%%10==0)
      cat(sprintf("Completed analysis of model %i of %i \n", i, N))
  }
  
  # Object manipulation
  lmat <- do.call(rbind, L1)
  tmat <- do.call(rbind,thoughtct)
  wmat <- do.call(rbind, wordct)
  wmod <- rowSums(wmat)
  tmod <- rowSums(tmat)
  lmod <- rowSums(lmat)
  
  # Taking away null values
  index <- which(!unlist(lapply(L1, is.null)))
  lb <- lb[index]
  
  # Retrieve semantic coherence scores
  semcoh <- c()
  for (z in seq(N))
    semcoh <- c(semcoh, mod.out$semcoh[z])
  semcoh <- semcoh[index]
  
  # Constructing output
  out <- list(N=N, K=K, glob.max=glob.max, lb=lb, lmat=lmat, tmat=tmat, wmat=wmat,
              lmod=lmod, tmod=tmod, wmod=wmod, semcoh=semcoh, L1mat=NULL, L1mod=NULL,
              mass.threshold=mass.threshold, cov.effects=NULL,
              var.matrix=NULL, confidence.ratings=NULL,
              align.global=align.global, reg.formula=reg.formula, 
              reg.nsims=reg.nsims, reg.parameter.index=reg.parameter.index)
  
  # Top of mass analysis
  if (mass.threshold != 1){
    
    # Initiate object
    bulk <- list()
    L1limited <- list()
    
    for(k in 1:K) {
      temp <- sort(ref.beta[k,], decreasing=TRUE)
      cum <- cumsum(temp)
      n <- min(which(cum>mass.threshold))
      bulk[[k]] <- order(ref.beta[k,], decreasing=TRUE)[1:n]
    }
    
    for(i in seq(N)) {
      
      # Distance in topic-word distributions
      betamod <- exp(mod.out$runout[[i]]$beta$logbeta[[1]])
      if(ncol(betamod)!=ncol(ref.beta)) next
      
      if(align.global==TRUE)
        align <- clue::solve_LSAP(ref.beta%*%t(betamod), maximum=TRUE)
      
      if(align.global==FALSE){
        product <- ref.beta%*%t(betamod)
        align <- apply(product,1,which.max)
      }
      
      betamod <- betamod[align,]
      
      lim <- vector(length=nrow(betamod))
      for(k in 1:nrow(betamod)) {
        lim[k] <- sum(abs(betamod[k,bulk[[k]]] - ref.beta[k, bulk[[k]]])) 
      }
      L1limited[[i]] <- lim
      
      if(verbose==TRUE & i%%10==0){
        cat(sprintf("Limited Mass Routine: Completed analysis of model %i of %i \n",
                    i, N))
      }
    }
    
    L1mat <- do.call(rbind, L1limited)
    L1mod <- rowMeans(L1mat)
    out$L1mat <- L1mat
    out$L1mod <- L1mod
  }
  
  # Confidence intervals abidance score for covariate effects
  if (!is.null(reg.formula)){
    var.matrix <- na.omit(var.matrix)
    confidence.ratings <- c()
    length(confidence.ratings) <- nrow(var.matrix)
    
    for(i in 1:nrow(var.matrix)){
      
      indicators <- c()
      length(indicators) <- K
      
      for(j in 1:K){
        if(var.matrix[i,j] > confidence.intervals[1,j] & 
             var.matrix[i,j] < confidence.intervals[2,j])
          indicators[j] <- 1
        else
          indicators[j] <- 0
      }
      
      score <- sum(indicators) / K
      confidence.ratings[i] <- score
    }
    
    out$cov.effects = cov.effects
    out$var.matrix = var.matrix
    out$confidence.ratings = confidence.ratings
  }
  
  # Load S3 Class and return
  class(out) <- "MultimodDiagnostic"
  return(out)
}

########################################
# Methods for MultimodDiagnostic class #
########################################

# Method Definitions
l1_topwords.MultimodDiagnostic <- function(obj) {
  smoothScatter(obj$wmat, obj$lmat, ylab="L1 Distance", 
                xlab="Top 10 Words in Common", main="L1 and Top Words")
}

l1_topdocs.MultimodDiagnostic <- function(obj) {
  smoothScatter(obj$tmat, obj$lmat, ylab="L1 Distance", 
                xlab="Top 10 Docs in Common", main="L1 and Top Documents")
}

topwords_topdocs.MultimodDiagnostic <- function(obj) {
  smoothScatter(obj$tmat, obj$wmat, xlab="Top 10 Docs in Common", 
                ylab="Top 10 Words in Common", main="Words and Documents")
}

bounds_words.MultimodDiagnostic <- function(obj) {
  plot(obj$lb, obj$wmod/1000, main="Bound and Words", 
       xlab="Bound", ylab="Proportion of 10 Words in Common")
  abline(lm(obj$wmod/1000 ~ obj$lb), lwd=2, col="blue")
  abline(h=.5, lty=2)
}

bounds_docs.MultimodDiagnostic <- function(obj) {
  plot(obj$lb, obj$tmod/1000, main="Bound and Docs", 
       xlab="Bound", ylab="Proportion of 10 Docs in Common")
  abline(lm(obj$tmod/1000 ~ obj$lb), lwd=2, col="blue")
  abline(h=.5, lty=2)
}

bounds_l1.MultimodDiagnostic <- function(obj) {
  plot(obj$lb, obj$lmod/100, main="Bound and L1", 
       xlab="Bound", ylab="Expected L1 Dist")
  abline(lm(obj$lmod/100 ~ obj$lb), lwd=2, col="blue")
  abline(h=1, lty=2)
}

wordsexpectation.MultimodDiagnostic <- function(obj) {
  hist(colMeans(obj$wmat), breaks=10, 
       xlab="Expected Common Words", main="Word Expectations over Runs")
}

docsexpectation.MultimodDiagnostic <- function(obj) {
  hist(colMeans(obj$tmat), breaks=10, 
       xlab="Expected Common Docs", main="Document Expectations over Runs")
}

l1expectation.MultimodDiagnostic <- function(obj) {
  hist(colMeans(obj$lmat), breaks=10, 
       xlab="Expected L1 Distance", main="L1 Expectations over Runs")
}

semcoh_words.MultimodDiagnostic <- function(obj) {
  y <- colMeans(obj$wmat)
  plot(obj$semcoh[[obj$glob.max]], y, xlab="Semantic Coherence", 
       ylab="Expected Words Shared", 
       main="Semantic Coherence and Words", ylim=c(0,10))
  abline(lm(y~obj$semcoh[[obj$glob.max]]), lwd=2, col="blue")
}

semcoh_docs.MultimodDiagnostic <- function(obj) {
  y <- colMeans(obj$tmat)
  plot(obj$semcoh[[obj$glob.max]], y, xlab="Semantic Coherence", 
       ylab="Expected Documents Shared",
       main="Semantic Coherence and Documents", ylim=c(0,10))
  abline(lm(y~obj$semcoh[[obj$glob.max]]), lwd=2, col="blue")
}

semcoh_l1.MultimodDiagnostic <- function(obj) {
  y <- colMeans(obj$lmat)
  plot(obj$semcoh[[obj$glob.max]], y, xlab="Semantic Coherence", 
       ylab="Expected L1 Distance",
       main="Semantic Coherence and L1 Distance", ylim=c(0,2))
  abline(lm(y~obj$semcoh[[obj$glob.max]]), lwd=2, col="blue")
}

l1_l1lim.MultimodDiagnostic <- function(obj) {
  smoothScatter(obj$lmat, obj$L1mat, ylab=paste("L1 Limited (", obj$mass.threshold, " mass)", sep=""), 
                xlab="L1", main="L1 and Limited L1")
}

topdocs_l1lim.MultimodDiagnostic <- function(obj) {
  smoothScatter(obj$wmat, obj$L1mat, ylab=paste("L1 Limited (", obj$mass.threshold, " mass)", sep=""), 
                xlab="Top 10 Documents in Common", main="Limited L1 and Top Documents")
}

bounds_l1lim.MultimodDiagnostic <- function(obj) {
  plot(obj$lb, obj$L1mod, main="Bound and Limited L1 Distance", xlab="Bound", 
       ylab=paste("Expected Limited L1 (", obj$mass.threshold, " mass) Distance", sep=""))
  abline(lm(obj$L1mod ~ obj$lb), lwd=2, col="blue")
  abline(h=.5, lty=2)
}

hist_l1lim.MultimodDiagnostic <- function(obj) {
  hist(colMeans(obj$L1mat), breaks=10, main="Histogram of Limited L1 by Topic", 
       xlab=paste("L1 Limited (", obj$mass.threshold, " mass)", sep=""))
}

semcoh_l1lim.MultimodDiagnostic <- function(obj) {
  y <- colMeans(obj$L1mat)
  plot(obj$semcoh[[obj$glob.max]], y, xlab="Semantic Coherence", 
       ylab=paste("Expected Limited L1 (", obj$mass.threshold, " mass) Distance", sep=""),
       main="Semantic Coherence and Limited L1 Distance", ylim=c(0,2))
  abline(lm(y~obj$semcoh[[obj$glob.max]]), lwd=2, col="blue")
}

topic_hist.MultimodDiagnostic <- function(obj, topics){
  oldpar <- par(no.readonly = T)
  N <- length(topics)
  root <- ceiling(sqrt(N))
  par(mfrow=c(root,root), oma=c(0,.5,1.5,0), mar=c(2,2,4,1))
  for(j in topics){
    hist(obj$var.matrix[, j], 
         main=paste("Topic", j),
         xlab="Coefficient Estimate",
         breaks=25)
    abline(v=0, lty=2)
    abline(v=obj$cov.effects[[obj$glob.max]]$parameters[[j]][[1]]$est[obj$reg.parameter.index],
           col=2)
  }
  title("Posterior Distribution of Covariate Effects By Topic", outer=TRUE)
  par(oldpar)
}

confidence_intervals.MultimodDiagnostic <- function(obj){
  hist(obj$confidence.ratings, 
       main="Distribution of .95 confidence-interval coverage for
       regression estimates",
       xlab="Confidence Interval Coverage")
}

plot.MultimodDiagnostic <- function(x, ind=NULL, topics=NULL, ...){
  obj <- x
  if(is.null(topics)) topics <- 1:obj$K
  if(is.null(ind)) ind <- 1:18
  .pardefault <- par(no.readonly = T)
  if(1 %in% ind) wordsexpectation.MultimodDiagnostic(obj, ...)
  if (length(ind) > 1) par(ask=T)
  if(2 %in% ind) docsexpectation.MultimodDiagnostic(obj, ...)
  if(!is.null(obj$cov.effect)){
    if(3 %in% ind) confidence_intervals.MultimodDiagnostic(obj, ...)
    if(4 %in% ind) topic_hist.MultimodDiagnostic(obj, topics, ...)
  }
  if(5 %in% ind) l1expectation.MultimodDiagnostic(obj, ...)
  if(6 %in% ind) l1_topwords.MultimodDiagnostic(obj, ...)
  if(7 %in% ind) l1_topdocs.MultimodDiagnostic(obj, ...)
  if(8 %in% ind) topwords_topdocs.MultimodDiagnostic(obj, ...)
  if(9 %in% ind) bounds_words.MultimodDiagnostic(obj, ...)
  if(10 %in% ind) bounds_docs.MultimodDiagnostic(obj, ...)
  if(11 %in% ind) bounds_l1.MultimodDiagnostic(obj, ...)
  if(12 %in% ind) semcoh_docs.MultimodDiagnostic(obj, ...)
  if(13 %in% ind) semcoh_l1.MultimodDiagnostic(obj, ...)
  if(14 %in% ind) semcoh_words.MultimodDiagnostic(obj, ...)
  if(!is.null(obj$L1mat)){
    if(15 %in% ind) hist_l1lim.MultimodDiagnostic(obj, ...)
    if(16 %in% ind) bounds_l1lim.MultimodDiagnostic(obj, ...)
    if(17 %in% ind) topdocs_l1lim.MultimodDiagnostic(obj, ...)
    if(18 %in% ind) semcoh_l1lim.MultimodDiagnostic(obj, ...)
  }
  par(.pardefault)
}

print.MultimodDiagnostic <- function(x,...){
  cat(paste("A multimodality diagnostic object for a collection of STM models.
            The diagnostic tests were run on", x$N, "models with", x$K, "topics."))
}

