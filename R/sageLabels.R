#' Displays verbose labels that describe topics and topic-covariate groups in
#' depth.
#' 
#' For each topic or, when there is a covariate at the bottom of the model, for
#' each topic-covariate group, sageLabels provides a list of the highest
#' marginal probability words, the highest marginal FREX words, the highest
#' marginal lift words, and the highest marginal score words, where marginal
#' means it is summing over all potential covariates.  It also provides each
#' topic's Kappa (words associated with each topic) and baselined Kappa
#' (baseline word distribution).
#' 
#' This can be used as an more detailed alternative to labelTopics.
#' 
#' @aliases sageLabels print.sageLabels
#' @param model A fitted STM model object.
#' @param n The number of words to print per topic/topic-covariate set. Default
#' is 7.
#' @return \item{marginal}{ A list of matrices, containing the high-probability
#' labels, FREX labels, lift labels, and high scoring words.  } \item{K}{ The
#' number of topics in the STM.  } \item{covnames}{ Names of the covariate
#' values used in the STM.  } \item{kappa}{Words associated with topics,
#' covariates, and topic/covariate interactions.} \item{kappa.m}{Baseline word
#' distribution.} \item{n}{ The n parameter passed by the user to this
#' function; number of words per topic or topic-covariate pair (when covariates
#' are used on the bottom of the model) } \item{cov.betas}{ Covariate-specific
#' beta matrices, listing for each covariate a matrix of highest-probability,
#' FREX, lift, and high scoring words.  Note that the actual vocabulary has
#' been substituted for word indices.  }
#' @export
sageLabels <- function(model, n=7) {
  logbeta <- model$beta$logbeta
  K <- model$settings$dim$K
  vocab <- model$vocab
  
  #Let's start by marginalizing
  margbeta <- exp(logbeta[[1]])
  if(length(logbeta) > 1) {
    weights <- model$settings$covariates$betaindex
    tab <- table(weights)
    weights <- tab/sum(tab)
    #marginalize
    margbeta <- margbeta*weights[1]
    for(i in 2:length(model$beta$logbeta)) {
      margbeta <- margbeta + exp(model$beta$logbeta[[i]])*weights[i]
    }
  }
  margbeta <- log(margbeta)
  #calculate all the usual suspects for the marginalized topics
  wordcounts <- model$settings$dim$wcounts$x
  
  frexlabels <- calcfrex(margbeta, .5, wordcounts)
  liftlabels <- calclift(margbeta, wordcounts)
  scorelabels <- calcscore(margbeta)
  problabels <- apply(margbeta, 1, order, decreasing=TRUE)
  
  
  frexlabels <- vocab[frexlabels[1:n,]]
  liftlabels <- vocab[liftlabels[1:n,]]
  scorelabels <- vocab[scorelabels[1:n,]]
  problabels <- vocab[problabels[1:n,]]
  
  # kappa.t
  labs <- lapply(model$beta$kappa$params, function(x) {
    windex <- order(x,decreasing=TRUE)[1:n]
    #we want to threshold by some minimal value.
    ifelse(x[windex]>1e-3, vocab[windex], "")
  }) 
  labs <- do.call(rbind, labs)
  labs <- labs[1:K,]
  
  # m + kappa.t
  labsm <- lapply(model$beta$kappa$params, function(x) {
    windex <- order(model$beta$kappa$m + x,decreasing=TRUE)[1:n]
    #we want to threshold by some minimal value.
    ifelse(x[windex]>1e-3, vocab[windex], "")
  }) 
  labsm <- do.call(rbind, labsm)
  labsm <- labsm[1:K,]
  
  # individual betas
  A <- model$settings$dim$A
  abetas <- vector(mode="list", length=A)
  for(a in 1:A) {
    dist <- model$beta$logbeta[[a]]
    abetas[[a]]$problabels <- matrix(vocab[apply(dist, 1, order, decreasing=TRUE)[1:n,]], nrow=K, ncol=n, byrow=TRUE)
    abetas[[a]]$frexlabels <- matrix(vocab[calcfrex(dist, .5, wordcounts)[1:n,]], nrow=K, ncol=n, byrow=TRUE)
    abetas[[a]]$liftlabels <- matrix(vocab[calclift(dist, wordcounts)[1:n,]], nrow=K, ncol=n, byrow=TRUE)
    abetas[[a]]$scorelabels <- matrix(vocab[calcscore(dist)[1:n,]], nrow=K, ncol=n, byrow=TRUE)
  }
    
  out <- list(marginal=list(frex=matrix(frexlabels, nrow=K, ncol=n, byrow=TRUE),
                            lift=matrix(liftlabels, nrow=K, ncol=n, byrow=TRUE),
                            score=matrix(scorelabels, nrow=K, ncol=n, byrow=TRUE), 
                            prob=matrix(problabels, nrow=K, ncol=n, byrow=TRUE)),
              kappa=labs,
              kappa.m=labsm,
              cov.betas=abetas, K=K, n=n, covnames=model$settings$covariates$yvarlevels)
  class(out) <- "sageLabels"
  return(out)
}

#' @method print sageLabels
#' @export
print.sageLabels <- function(x, ...) {
  topicnums <- 1:x$K
  
  #copying old stuff below
  for(i in topicnums) {
    toprint <- sprintf("Topic %i: \n \t Marginal Highest Prob: %s \n \t Marginal FREX: %s \n \t Marginal Lift: %s \n \t Marginal Score: %s \n \n", 
                       i, 
                       commas(x$marginal$prob[i,]),
                       commas(x$marginal$frex[i,]),
                       commas(x$marginal$lift[i,]),
                       commas(x$marginal$score[i,]))
    toprint2 <- sprintf(" \t Topic Kappa: %s \n \t Kappa with Baseline: %s \n \n",
                        commas(x$kappa[i,]),
                        commas(x$kappa.m[i,]))
    combine <- paste0(toprint, toprint2)
    for(a in 1:length(x$cov.betas)) {
      text <- sprintf(" \t Covariate %s: \n \t \t Marginal Highest Prob: %s \n \t \t Marginal FREX: %s \n \t \t Marginal Lift: %s \n \t \t Marginal Score: %s \n", 
                      x$covnames[a], 
                      commas(x$cov.betas[[a]]$problabels[i,]),
                      commas(x$cov.betas[[a]]$frexlabels[i,]),
                      commas(x$cov.betas[[a]]$liftlabels[i,]),
                      commas(x$cov.betas[[a]]$scorelabels[i,]))
      combine <- paste0(combine,text)
    }
    cat(combine)
  }
}