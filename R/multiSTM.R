#########################################################################
#         Generalized Model Stability Analysis Routines for STM         #
#             Antonio Coppola, Harvard University, July 2014            #
#########################################################################

#' Analyze Stability of Local STM Mode
#' 
#' This function performs a suite of tests aimed at assessing the global
#' behavior of an STM model, which may have multiple modes. The function takes
#' in a collection of differently initialized STM fitted objects and selects a
#' reference model against which all others are benchmarked for stability. The
#' function returns an output of S3 class 'MultimodDiagnostic', with associated
#' plotting methods for quick inspection of the test results.
#' 
#' The purpose of this function is to automate and generalize the stability
#' analysis routines for topic models that are introduced in Roberts, Margaret
#' E., Brandon M. Stewart, and Dustin Tingley: "Navigating the Local Modes of
#' Big Data: The Case of Topic Models" (2014). For more detailed discussion
#' regarding the background and motivation for multimodality analysis, please
#' refer to the original article. See also the documentation for
#' \code{\link{plot.MultimodDiagnostic}} for help with the plotting methods
#' associated with this function.
#' 
#' @aliases multiSTM print.MultimodDiagnostic
#' @param mod.out The output of a \code{selectModel()} run. This is a list of
#' model outputs the user has to choose from, which all take the same form as
#' the output from a STM model. Currently only works with models without
#' content covariates.
#' @param ref.model An integer referencing the element of the list in
#' \code{mod.out} which contains the desired reference model.  When set to the
#' default value of \code{NULL} this chooses the model with the largest value
#' of the approximate variational bound.
#' @param align.global A boolean parameter specifying how to align the topics
#' of two different STM fitted models. The alignment is performed by solving
#' the linear sum assignment problem using the Hungarian algorithm. If
#' \code{align.global} is set to \code{TRUE}, the Hungarian algorithm is run
#' globally on the topic-word matrices of the two models that are being
#' compared. The rows of the matrices are aligned such as to minimize the sum
#' of their inner products. This results in each topic in the current runout
#' being matched to a unique topic in the reference model. If
#' \code{align.global} is, conversely, set to \code{FALSE}, the alignment
#' problem is solved locally. Each topic in the current runout is matched to
#' the one topic in the reference models that yields minimum inner product.
#' This means that multiple topics in the current runout can be matched to a
#' single topic in the reference model, and does not guarantee that all the
#' topics in the reference model will be matched.
#' @param mass.threshold A parameter specifying the portion of the probability
#' mass of topics to be used for model analysis. The tail of the probability
#' mass is disregarded accordingly. If \code{mass.threshold} is different from
#' 1, both the full-mass and partial-mass analyses are carried out.
#' @param reg.formula A formula for estimating a regression for each model in
#' the ensemble, where the documents are the units, the outcome is the
#' proportion of each document about a topic in an STM model, and the
#' covariates are the document-level metadata. The formula should have an
#' integer or a vector of numbers on the left-hand side, and an equation with
#' covariates on the right-hand side. If the left-hand side is left blank, the
#' regression is performed on all topics in the model. The formula is
#' exclusively used for building calls to \code{estimateEffect()}, so see the
#' documentation for \code{estimateEffect()} for greater detail about the
#' regression procedure. If \code{reg.formula} is null, the covariate effect
#' stability analysis routines are not performed. The regressions incorporate
#' uncertainty by using an approximation to the average covariance matrix
#' formed using the global parameters.
#' @param metadata A dataframe where the predictor variables in
#' \code{reg.formula} can be found. It is necessary to include this argument if
#' \code{reg.formula} is specified.
#' @param reg.nsims The number of simulated draws from the variational
#' posterior for each call of \code{estimateEffect()}. Defaults to 100.
#' @param reg.parameter.index If \code{reg.formula} is specified, the function
#' analyzes the stability across runs of the regression coefficient for one
#' particular predictor variable. This argument specifies which predictor
#' variable is to be analyzed. A value of 1 corresponds to the intercept, a
#' value of 2 correspond to the first predictor variable in \code{reg.formula},
#' and so on. Support for multiple concurrent covariate effect stability
#' analyses is forthcoming.
#' @param verbose If set to \code{TRUE}, the function will report progress.
#' @param from.disk If set to \code{TRUE}, \code{multiSTM()} will load the
#' input models from disk rather than from RAM. This option is particularly
#' useful for dealing with large numbers of models, and is intended to be used
#' in conjunction with the \code{to.disk} option of \code{selectModel()}.
#' \code{multiSTM()} inspects the current directory for RData files.
#' @return An object of 'MultimodDiagnostic' S3 class, consisting of a list
#' with the following components: \item{N}{The number of fitted models in the
#' list of model outputs that was supplied to the function for the purpose of
#' stability analysis.} \item{K}{The number of topics in the models.}
#' \item{glob.max}{The index of the reference model in the list of model
#' outputs (\code{mod.out}) that was supplied to the function. The reference
#' model is selected as the one with the maximum bound value at convergence.}
#' \item{lb}{A list of the maximum bound value at convergence for each of the
#' fitted models in the list of model outputs. The list has length N.}
#' \item{lmat}{A K-by-N matrix reporting the L1-distance of each topic from the
#' corresponding one in the reference model. This is defined as:
#' \deqn{L_{1}=\sum_{v}|\beta_{k,v}^{ref}-\beta_{k,v}^{cand}|} Where the beta
#' matrices are the topic-word matrices for the reference and the candidate
#' model.} \item{tmat}{A K-by-N matrix reporting the number of "top documents"
#' shared by the reference model and the candidate model. The "top documents"
#' for a given topic are defined as the 10 documents in the reference corpus
#' with highest topical frequency.} \item{wmat}{A K-by-N matrix reporting the
#' number of "top words" shared by the reference model and the candidate model.
#' The "top words" for a given topic are defined as the 10 highest-frequency
#' words.} \item{lmod}{A vector of length N consisting of the row sums of the
#' \code{lmat} matrix.} \item{tmod}{A vector of length N consisting of the row
#' sums of the \code{tmat} matrix.} \item{wmod}{A vector of length N consisting
#' of the row sums of the \code{wmat} matrix.} \item{semcoh}{Semantic coherence
#' values for each topic within each model in the list of model outputs.}
#' \item{L1mat}{A K-by-N matrix reporting the limited-mass L1-distance of each
#' topic from the corresponding one in the reference model. Similar to
#' \code{lmat}, but computed using only the top portion of the probability mass
#' for each topic, as specified by the \code{mass.threshol} parameter.
#' \code{NULL} if \code{mass.treshold==1}.} \item{L1mod}{A vector of length N
#' consisting of the row means of the \code{L1mat} matrix.}
#' \item{mass.threshold}{The mass threshold argument that was supplied to the
#' function.} \item{cov.effects}{A list of length N containing the output of
#' the run of \code{estimateEffect()} on each candidate model with the given
#' regression formula. \code{NULL} if no regression formula is given.}
#' \item{var.matrix}{A K-by-N matrix containing the estimated variance for each
#' of the fitted regression parameters. \code{NULL} if no regression formula is
#' given.} \item{confidence.ratings}{A vector of length N, where each entry
#' specifies the proportion of regression coefficient estimates in a candidate
#' model that fall within the .95 confidence interval for the corresponding
#' estimate in the reference model.} \item{align.global}{The alignment control
#' argument that was supplied to the function.} \item{reg.formula}{The
#' regression formula that was supplied to the function.} \item{reg.nsims}{The
#' \code{reg.nsims} argument that was supplied to the function.}
#' \item{reg.parameter.index}{The \code{reg.parameter.index} argument that was
#' supplied to the function.}
#' @author Antonio Coppola (Harvard University), Brandon Stewart (Princeton
#' University), Dustin Tingley (Harvard University)
#' @seealso \code{\link{plot.MultimodDiagnostic}} \code{\link{selectModel}}
#' \code{\link{estimateEffect}}
#' @references Roberts, M., Stewart, B., & Tingley, D. (2016).
#' "Navigating the Local Modes of Big Data: The Case of Topic Models. In Data
#' Analytics in Social Science, Government, and Industry." New York: Cambridge
#' University Press.
#' @keywords stm multimodality
#' @examples
#' 
#' \dontrun{
#' 
#' # Example using Gadarian data
#' temp<-textProcessor(documents=gadarian$open.ended.response, 
#'                     metadata=gadarian)
#' meta<-temp$meta
#' vocab<-temp$vocab
#' docs<-temp$documents
#' out <- prepDocuments(docs, vocab, meta)
#' docs<-out$documents
#' vocab<-out$vocab
#' meta <-out$meta
#' set.seed(02138)
#' mod.out <- selectModel(docs, vocab, K=3, 
#'                        prevalence=~treatment + s(pid_rep), 
#'                        data=meta, runs=20)
#' 
#' out <- multiSTM(mod.out, mass.threshold = .75, 
#'                 reg.formula = ~ treatment,
#'                 metadata = gadarian)
#' plot(out)
#' 
#' # Same example as above, but loading from disk
#' mod.out <- selectModel(docs, vocab, K=3, 
#'                        prevalence=~treatment + s(pid_rep), 
#'                        data=meta, runs=20, to.disk=T)
#' 
#' out <- multiSTM(from.disk=T, mass.threshold = .75, 
#'                 reg.formula = ~ treatment,
#'                 metadata = gadarian)
#' }
#' @export
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

#' Plotting Method for Multimodality Diagnostic Objects
#' 
#' The plotting method for objects of the S3 class 'MultimodDiagnostic', which
#' are returned by the function \code{multiSTM()}, which performs a battery of
#' tests aimed at assessing the stability of the local modes of an STM model.
#' 
#' This methods generates a series of plots, which are indexed as follows. If a
#' subset of the plots is required, specify their indexes using the \code{ind}
#' argument. Please note that not all plot types are available for every object
#' of class 'MultimodDiagnostic': \enumerate{ \item Histogram of Expected
#' Common Words: Generates a 10-bin histogram of the column means of
#' \code{obj$wmat}, a K-by-N matrix reporting the number of "top words" shared
#' by the reference model and the candidate model. The "top words" for a given
#' topic are defined as the 10 highest-frequency words.  \item Histogram of
#' Expected Common Documents: Generates a 10-bin histogram of the column means
#' of \code{obj$tmat}, a K-by-N matrix reporting the number of "top documents"
#' shared by the reference model and the candidate model. The "top documents"
#' for a given topic are defined as the 10 documents in the reference corpus
#' with highest topical frequency.  \item Distribution of .95
#' Confidence-Interval Coverage for Regression Estimates: Generates a histogram
#' of \code{obj$confidence.ratings}, a vector whose entries specify the
#' proportion of regression coefficient estimates in a candidate model that
#' fall within the .95 confidence interval for the corresponding estimate in
#' the reference model. This can only be generated if
#' \code{obj$confidence.ratings} is non-\code{NULL}.  \item Posterior
#' Distributions of Covariate Effect Estimates By Topic: Generates a square
#' matrix of plots, each depicting the posterior distribution of the regression
#' coefficients for the covariate specified in \code{obj$reg.parameter.index}
#' for one topic. The topics for which the plots are to be generated are
#' specified by the \code{topics} argument. If the length of \code{topics} is
#' not a perfect square, the plots matrix will include white space. The plots
#' have a dashed black vertical line at zero, and a continuous red vertical
#' line indicating the coefficient estimate in the reference model. This can
#' only be generated if \code{obj$cov.effects} is non-\code{NULL}.  \item
#' Histogram of Expected L1-Distance From Reference Model: Generates a 10-bin
#' histogram of the column means of \code{obj$lmat}, a K-by-N matrix reporting
#' the L1-distance of each topic from the corresponding one in the reference
#' model.  \item L1-distance vs. Top-10 Word Metric: Produces a smoothed color
#' density representation of the scatterplot of \code{obj$lmat} and
#' \code{obj$wmat}, the metrics for L1-distance and shared top-words, obtained
#' through a kernel density estimate. This can be used to validate the metrics
#' under consideration.  \item L1-distance vs. Top-10 Docs Metric: Produces a
#' smoothed color density representation of the scatterplot of \code{obj$lmat}
#' and \code{obj$tmat}, the metrics for L1-distance and shared top-documents,
#' obtained through a kernel density estimate. This can be used to validate the
#' metrics under consideration.  \item Top-10 Words vs. Top-10 Docs Metric:
#' Produces a smoothed color density representation of the scatterplot of
#' \code{obj$wmat} and \code{obj$tmat}, the metrics for shared top-words and
#' shared top-documents, obtained through a kernel density estimate. This can
#' be used to validate the metrics under consideration.  \item Maximized Bound
#' vs. Aggregate Top-10 Words Metric: Generates a scatter plot with linear
#' trendline for the maximized bound vector (\code{obj$lb}) and a linear
#' transformation of the top-words metric aggregated by model
#' (\code{obj$wmod/1000}).  \item Maximized Bound vs. Aggregate Top-10 Docs
#' Metric: Generates a scatter plot with linear trendline for the maximized
#' bound vector (\code{obj$lb}) and a linear transformation of the top-docs
#' metric aggregated by model (\code{obj$tmod/1000}).  \item Maximized Bound
#' vs. Aggregate L1-Distance Metric: Generates a scatter plot with linear
#' trendline for the maximized bound vector (\code{obj$lb}) and a linear
#' transformation of the L1-distance metric aggregated by model
#' (\code{obj$tmod/1000}).  \item Top-10 Docs Metric vs. Semantic Coherence:
#' Generates a scatter plot with linear trendline for the reference-model
#' semantic coherence scores and the column means of \code{object$tmat}.  \item
#' L1-Distance Metric vs. Semantic Coherence: Generates a scatter plot with
#' linear trendline for the reference-model semantic coherence scores and the
#' column means of \code{object$lmat}.  \item Top-10 Words Metric vs. Semantic
#' Coherence: Generates a scatter plot with linear trendline for the
#' reference-model semantic coherence scores and the column means of
#' \code{object$wmat}.  \item Same as \code{5}, but using the limited-mass
#' L1-distance metric. Can only be generated if \code{obj$mass.threshold != 1}.
#' \item Same as \code{11}, but using the limited-mass L1-distance metric. Can
#' only be generated if \code{obj$mass.threshold != 1}.  \item Same as
#' \code{7}, but using the limited-mass L1-distance metric. Can only be
#' generated if \code{obj$mass.threshold != 1}.  \item Same as \code{13}, but
#' using the limited-mass L1-distance metric. Can only be generated if
#' \code{obj$mass.threshold != 1}.  }
#' 
#' @param x An object of S3 class 'MultimodDiagnostic'. See
#' \code{\link{multiSTM}}.
#' @param ind An integer of list of integers specifying which plots to generate
#' (see details). If \code{NULL} (default), all plots are generated.
#' @param topics An integer or vector of integers specifying the topics for
#' which to plot the posterior distribution of covariate effect estimates. If
#' \code{NULL} (default), plots are generated for every topic in the S3 object.
#' @param ... Other arguments to be passed to the plotting functions.
#' @author Brandon M. Stewart (Princeton University) and Antonio Coppola
#' (Harvard University)
#' @seealso \code{\link{multiSTM}}
#' @references Roberts, M., Stewart, B., & Tingley, D. (Forthcoming).
#' "Navigating the Local Modes of Big Data: The Case of Topic Models. In Data
#' Analytics in Social Science, Government, and Industry." New York: Cambridge
#' University Press.
#' @keywords stm multimodality
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' # Example using Gadarian data
#' 
#' temp<-textProcessor(documents=gadarian$open.ended.response, 
#'                     metadata=gadarian)
#' meta<-temp$meta
#' vocab<-temp$vocab
#' docs<-temp$documents
#' out <- prepDocuments(docs, vocab, meta)
#' docs<-out$documents
#' vocab<-out$vocab
#' meta <-out$meta
#' set.seed(02138)
#' mod.out <- selectModel(docs, vocab, K=3, 
#'                        prevalence=~treatment + s(pid_rep), 
#'                        data=meta, runs=20)
#' 
#' out <- multiSTM(mod.out, mass.threshold = .75, 
#'                 reg.formula = ~ treatment,
#'                 metadata = gadarian)
#' 
#' plot(out)
#' plot(out, 1)
#' }
#' @export
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

#' @export
#' @method print MultimodDiagnostic
print.MultimodDiagnostic <- function(x,...){
  cat(paste("A multimodality diagnostic object for a collection of STM models.
            The diagnostic tests were run on", x$N, "models with", x$K, "topics."))
}