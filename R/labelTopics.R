#Topic Labeling according to a series of metrics.

labelTopics <- function (model, topics=NULL, n = 7, frexweight=.5) {
  logbeta <- model$beta$logbeta
  K <- model$settings$dim$K
  vocab <- model$vocab
  if(is.null(topics)) topics <- 1:nrow(logbeta[[1]])    
  #make a switch for presence of contnet covariate
  aspect <- length(logbeta)>1

  out <- list()
  if(!aspect) {    
    logbeta <- logbeta[[1]]
    wordcounts <- model$settings$dim$wcounts$x
    #Calculate FREX Score
    frexlabels <- calcfrex(logbeta, frexweight, wordcounts)
    
    #Calculate Lift (Taddys thing this is beta_k,v divided by the empirical term probability)
    liftlabels <- calclift(logbeta, wordcounts)

    #Calculate score (Chang in LDA package etc.)
    scorelabels <- calcscore(logbeta)
    
    #standard labels
    problabels <- apply(logbeta, 1, order, decreasing=TRUE)
    
    for(k in 1:K) {
      out$prob[[k]] <- vocab[problabels[1:n,k]]
      out$frex[[k]] <- vocab[frexlabels[1:n,k]]
      out$lift[[k]] <- vocab[liftlabels[1:n,k]]
      out$score[[k]] <- vocab[scorelabels[1:n,k]]
    }
    out <- lapply(out, do.call, what=rbind)
  } else {
    labs <- lapply(model$beta$kappa$params, function(x) {
      windex <- order(x,decreasing=TRUE)[1:n]
      ifelse(x[windex]>0, vocab[windex], "")
    }) 
    labs <- do.call(rbind, labs)
    A <- model$settings$dim$A
    anames <- model$settings$covariates$yvarlevels
    i1 <- K + 1; i2 <- K + A; 
    
    out$topics <- labs[1:K,]
    out$covariate <- labs[i1:i2,]
    rownames(out$covariate) <- anames
    out$interaction <- labs[(i2+1):nrow(labs),]
  }
  out$topicnums <- topics
  class(out) <- "labelTopics"
  return(out)
}

print.labelTopics <- function(x,...) {  
  #test of its an aspect model or not
  if(names(x)[1]!="topics") {
    #Standard Non-aspect models
    K <- nrow(x$prob)  
    for(k in x$topicnums) {
      toprint <- sprintf("Topic %i Top Words:\n \t Highest Prob: %s \n \t FREX: %s \n \t Lift: %s \n \t Score: %s \n", 
                         k, 
                         commas(x$prob[k,]),
                         commas(x$frex[k,]),
                         commas(x$lift[k,]),
                         commas(x$score[k,]))
      cat(toprint)
    }
  } else {
    K <- nrow(x$topics)
    covlevels <- rownames(x$covariate)
    #Aspect Models 
    topiclabs <- c("Topic Words:\n",sprintf("Topic %i: %s \n", 
                                            x$topicnums, apply(x$topics[x$topicnums,], 1, commas)))  
    aspects <- c("Covariate Words:\n", sprintf("Group %s: %s \n", 
                                               covlevels, 
                                               apply(x$covariate, 1, commas)))
    #interactions are labeled along a sequence of 1:K, 1:K etc.
    interactions <- c("Topic-Covariate Interactions:\n") 
    intlabs <- apply(x$interaction,1,commas)
    for(k in 1:x$topicnums) {
      wseq <- seq(from=k, by=K, to=length(intlabs))
      out <- sprintf("Topic %i, Group %s: %s \n", k, covlevels, intlabs[wseq])
      interactions <- c(interactions,out,"\n")
    }
    labels <- c(topiclabs,"\n",aspects,"\n",interactions)
    cat(labels)
  }
}
  
  
