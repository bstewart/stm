#' Label topics
#' 
#' Generate a set of words describing each topic from a fitted STM object.
#' Uses a variety of labeling algorithms (see details).
#' 
#' Four different types of word weightings are printed with label topics.
#' 
#' Highest Prob: are the words within each topic with the highest probability
#' (inferred directly from topic-word distribution parameter \eqn{\beta}).
#' 
#' FREX: are the words that are both frequent and exclusive, identifying words
#' that distinguish topics.  This is calculated by taking the harmonic mean of
#' rank by probability within the topic (frequency) and rank by distribution of
#' topic given word \eqn{p(z|w=v)} (exclusivity).  In estimating exclusivity we
#' use a James-Stein type shrinkage estimator of the distribution
#' \eqn{p(z|w=v)}.  More information can be found in the documentation for the
#' internal function \code{\link{calcfrex}} and \code{\link{js.estimate}}.
#' 
#' Score and Lift are measures provided in two other popular text mining
#' packages. For more information on type Score, see the R package
#' \code{\link{lda}} or the internal function \code{\link{calcscore}}.  
#' For more information on type Lift, see the R package \code{maptpx}
#' or or the internal function \code{\link{calclift}}.
#' 
#' @aliases labelTopics print.labelTopics
#' @param model An \code{STM} model object.
#' @param topics A vector of numbers indicating the topics to include.  Default
#' is all topics.
#' @param n The desired number of words (per type) used to label each topic.
#' Must be 1 or greater.
#' @param frexweight A weight used in our approximate FREX scoring algorithm
#' (see details).
#' @return A labelTopics object (list) \item{prob }{matrix of highest
#' probability words} \item{frex }{matrix of highest ranking frex words}
#' \item{lift }{matrix of highest scoring words by lift} \item{score }{matrix
#' of best words by score} \item{topicnums }{a vector of topic numbers which
#' correspond to the rows}
#' @seealso \code{\link{stm}} \code{\link{plot.STM}} 
#' \code{\link{calcfrex}} \code{\link{js.estimate}} \code{\link{calcscore}} \code{\link{calclift}}
#' @examples
#' labelTopics(gadarianFit)
#' @export
labelTopics <- function (model, topics=NULL, n = 7, frexweight=.5) {
  if(n<1) stop("n must be 1 or greater")
  logbeta <- model$beta$logbeta
  K <- model$settings$dim$K
  vocab <- model$vocab
  if(is.null(topics)) topics <- 1:nrow(logbeta[[1]])    
  #make a switch for presence of contnet covariate
  aspect <- length(logbeta)>1

  out <- list()
  if(!aspect) {
    out$prob <- list()
    out$frex <- list()
    out$lift <- list()
    out$score <- list()
    logbeta <- logbeta[[1]]
    wordcounts <- model$settings$dim$wcounts$x
    #Calculate FREX Score
    frexlabels <- try(calcfrex(logbeta, frexweight, wordcounts),silent=TRUE)

    #Calculate Lift (Taddys thing this is beta_k,v divided by the empirical term probability)
    liftlabels <- try(calclift(logbeta, wordcounts), silent=TRUE)

    #Calculate score (Chang in LDA package etc.)
    scorelabels <- try(calcscore(logbeta), silent=TRUE)
    
    #standard labels
    problabels <- apply(logbeta, 1, order, decreasing=TRUE)
    
    for(k in 1:K) {
      out$prob[[k]] <- vocab[problabels[1:n,k]]
      if(inherits(frexlabels,"try-error")) {
        out$frex[[k]] <- "FREX encountered an error and failed to run"
      } else {
        out$frex[[k]] <- vocab[frexlabels[1:n,k]]        
      }
      if(inherits(liftlabels,"try-error")) {
        out$lift[[k]] <- "Lift encountered an error and failed to run"
      } else {
        out$lift[[k]] <- vocab[liftlabels[1:n,k]]        
      }  
      if(inherits(scorelabels,"try-error")) {
        out$lift[[k]] <- "Score encountered an error and failed to run"
      } else {
        out$score[[k]] <- vocab[scorelabels[1:n,k]]        
      }  
    }
    out <- lapply(out, do.call, what=rbind)
  } else {
    labs <- lapply(model$beta$kappa$params, function(x) {
      windex <- order(x,decreasing=TRUE)[1:n]
      #we want to threshold by some minimal value.
      ifelse(x[windex]>1e-3, vocab[windex], "")
    }) 
    labs <- do.call(rbind, labs)
    A <- model$settings$dim$A
    anames <- model$settings$covariates$yvarlevels
    i1 <- K + 1; i2 <- K + A; 
    intnums <- (i2+1):nrow(labs)
    
    out$topics <- labs[topics,,drop=FALSE]
    out$covariate <- labs[i1:i2,,drop=FALSE]
    rownames(out$covariate) <- anames
    if(model$settings$kappa$interactions) {
      tindx <- rep(1:K, each=A)
      intnums <- intnums[tindx %in% topics]
      out$interaction <- labs[intnums,,drop=FALSE]
    }
  }
  out$topicnums <- topics
  class(out) <- "labelTopics"
  return(out)
}

#' @method print labelTopics
#' @export
print.labelTopics <- function(x,...) {  
  #test of its an aspect model or not
  if(names(x)[1]!="topics") {
    #Standard Non-aspect models 
    for(i in 1:length(x$topicnums)) {
      toprint <- sprintf("Topic %i Top Words:\n \t Highest Prob: %s \n \t FREX: %s \n \t Lift: %s \n \t Score: %s \n", 
                         x$topicnums[i], 
                         commas(x$prob[x$topicnums[i],]),
                         commas(x$frex[x$topicnums[i],]),
                         commas(x$lift[x$topicnums[i],]),
                         commas(x$score[x$topicnums[i],]))
      cat(toprint)
    }
  } else {
    covlevels <- rownames(x$covariate)
    #Aspect Models 
    topiclabs <- c("Topic Words:\n",sprintf("Topic %i: %s \n", 
                                            x$topicnums, apply(x$topics, 1, commas)))  
    aspects <- c("Covariate Words:\n", sprintf("Group %s: %s \n", 
                                               covlevels, 
                                               apply(x$covariate, 1, commas)))
    if("interaction" %in% names(x)) {
      interactions <- c("Topic-Covariate Interactions:\n") 
      intlabs <- apply(x$interaction,1,commas)
      topicnums <- rep(x$topicnums, times=length(covlevels))
      
      for(i in x$topicnums) {
        out <- sprintf("Topic %i, Group %s: %s \n", i, covlevels, intlabs[which(topicnums==i)])
        interactions <- c(interactions,out,"\n")
      }
      labels <- c(topiclabs,"\n",aspects,"\n",interactions)
    } else {
      labels <- c(topiclabs, "\n", aspects)
    }
    cat(labels)
  }
}