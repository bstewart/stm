#Functions for plotting STM objects

plot.STM <- function(x, 
                     type=c("summary", "labels", "perspectives", "hist"),
                     n=NULL, topics=NULL,
                     labeltype=c("prob", "frex", "lift", "score"), 
                     frexw=.5,
                     main=NULL, xlim=NULL, ylim=NULL, xlab=NULL, family="",
                     width=80, 
                     covarlevels=NULL, plabels=NULL, text.cex=1, custom.labels=NULL,
                     topic.names=NULL,
                     ...){
  model <- x
  type <- match.arg(type)
  
  contentcov <- length(model$beta$logbeta)!=1
  if(contentcov & !missing(labeltype)) stop("Cannot specify label type for content covariate models.")
  
  labeltype <- match.arg(labeltype)
  if(!is.null(custom.labels)) labeltype <- "custom"
  if(is.null(n)) n <- switch(type, 
                             summary=3, 
                             labels=20,
                             perspectives=25,
                             hist=3)
  
  if(type!="perspectives" & is.null(topics)) topics <- 1:model$settings$dim$K
  
  if(labeltype!="custom"){
    if(type != "perspectives") {
      lab <- labelTopics(model, topics=topics, n = n, frexweight=frexw)
      if(contentcov) {
        lab <- lab$topics
      } else {
        lab <- lab[[labeltype]]
      }    
    }
  } else {
    lab <- custom.labels
    if(length(lab)!=length(topics)) lab <- rep_len(lab, length.out=length(topics))
  }
  
  if(!is.null(topic.names))  topic.names <- rep_len(topic.names, length.out=length(topics))
  
  ##############
  # Summary Plot
  ##############
  if(type=="summary") {
    if(labeltype!="custom"){
      lab <- apply(lab, 1, commas)
      if(!contentcov) lab <- lab[topics]
    }
    if(is.null(topic.names)) {
      topic.names <- sprintf("Topic %i:", topics)
    }   
    lab <- sprintf("%s %s", topic.names, lab)
    
    frequency <- colMeans(model$theta[,topics])
    invrank <- order(frequency, decreasing=FALSE)
    if(is.null(xlim)) xlim <- c(0,min(2*max(frequency), 1))
    if(is.null(ylim)) ylim <- c(0,length(topics))
    if(is.null(main)) main <- "Top Topics"
    if(is.null(xlab)) xlab <- "Expected Topic Proportions"
    ylab <- ""

    plot(c(0,0), type="n", xlim=xlim, ylim=ylim, main=main, 
         yaxt="n", 
         ylab=ylab, xlab=xlab, ...)
    for(i in 1:length(invrank)) {
      lines(c(0,frequency[invrank[i]]), c(i, i))
      text(frequency[invrank[i]]+.01, i , lab[invrank[i]],family=family,pos=4, cex=text.cex)
    }
  }
  
  ##############
  # Labels
  ##############
  if(type=="labels") {
    if(contentcov) {
      #we want to marginalize over the aspects here
      weights <- model$settings$covariates$betaindex
      tab <- table(weights)
      weights <- tab/sum(tab)
      #marginalize
      beta <- exp(model$beta$logbeta[[1]])*weights[1]
      for(i in 2:length(model$beta$logbeta)) {
        beta <- beta + exp(model$beta$logbeta[[i]])*weights[i]
      }
      lab <- t(apply(beta, 1, function(x) model$vocab[order(x, decreasing=TRUE)[1:n]]))
    }
    plot(c(0,0), type="n", main=main,
         ylim=c((length(topics)+.9),1), #note- this flips the coordinate system
         xlim=c(1,n+1), 
         xaxt="n", yaxt="n", xlab="", ylab="",...)
    
    if(labeltype!="custom") {
      lab <- apply(lab, 1, commas)
      lab <- lab[topics]
    }
    if(is.null(topic.names)) {
      topic.names <- sprintf("Topic %i:", topics)
    }
    lab <- lapply(lab, strwrap, width=width)

    
    for(i in 1:length(topics)){
      if(i!=length(topics)) lines(c(-1,n+3), c(i+1,i+1), lty=2)
      string <- paste0(lab[[i]],collapse="\n")
      string <- sprintf("%s \n %s", topic.names[i], string)
      text(n/2, (i + .5), string, family=family, cex=text.cex)
    }
  }
  
  ##############
  # Perspectives
  ##############
  if(type=="perspectives") {
    if(!contentcov) covarlevels <- c(1,1)
    
    #check the topics
    if(is.null(topics)) stop("Must specify one or two topic numbers using topics")
    if(length(topics)>2) stop("Too many topics specified.")
    if(length(topics)==1) {
      topics <- rep(topics,2)
    }
    sametopics <- topics[1]==topics[2]

    #check the aspects
    if(is.null(covarlevels)) {
      covarlevels <- c(1,2)
    } else {
      if(length(covarlevels)>2) stop("More than two covariate levels specified.")
      if(length(covarlevels)==1) covarlevels <- rep(covarlevels,2)
      covlabels <-  model$seetings$covariates$yvarlevels
      if(is.character(covarlevels)) {
        covarlevels <- pmatch(covarlevels, model$settings$covariates$yvarlevels)
        if(any(is.na(covarlevels))) stop("Unrecognized covariate levels")
      } 
    }
    samecovars <- covarlevels[1]==covarlevels[2]
    
    #determine what the axis labels should be
    if(is.null(plabels)) {
      plabels <- vector(length=2)
      if(!contentcov) {
        plabels[1] <- paste0(c("Topic ", topics[1]),collapse="")
        plabels[2] <- paste0(c("Topic ", topics[2]),collapse="")
      } else {
        if(sametopics & !samecovars) {
          plabels[1] <- paste0(c(model$settings$covariates$yvarlevels[covarlevels[1]], 
                              "\n", "(Topic ", topics[1], ")"),
                            collapse="")
          plabels[2] <- paste0(c(model$settings$covariates$yvarlevels[covarlevels[2]], 
                                 "\n", "(Topic ", topics[2], ")"),
                               collapse="")
        } else {
          plabels[1] <- paste0(c("Topic ", topics[1], ") \n",
                              "(", 
                              model$settings$covariates$yvarlevels[covarlevels[1]], ")"),
                            collapse="")
          plabels[2] <- paste0(c("Topic ", topics[2], ") \n",
                                 "(", 
                                 model$settings$covariates$yvarlevels[covarlevels[2]], ")"),
                               collapse="")
        }
      }
    }

    
    #Retrieve the relevant logbetas
    # okay next thing to try is regularized estimates of the ratio.
    left <- model$beta$logbeta[[covarlevels[1]]][topics[1],]
    right <- model$beta$logbeta[[covarlevels[2]]][topics[2],]
    nd2 <- floor(n/2)
    words <- unique(c(order(left, decreasing=TRUE)[1:nd2], order(right,decreasing=TRUE)[1:nd2]))
    
    #Take the probability for the topic they are the highest in.  
    # Scale it so the max is 4 times as large, and then take the inverse hyperbolic sine func.
    # this basically has the effect of flattening it out a bit so the giant things aren't huge.
    scale <- pmax(exp(left[words]) + exp(right[words]))
    scale <- asinh((scale/max(scale))*4)
    
    diff <- exp(left[words]) - exp(right[words])
    diff <- diff/max(abs(diff))
    #Create the plot
    plot(c(0,0), xlim=c(max(diff)+.1,min(diff) - .1), ylim=c(-2*length(diff),6*length(diff)), type="n", 
         xaxt="n", xlab="", main=main, yaxt="n", ylab="", bty="n",...)
    segments(0,0, 0, -2*length(diff),lty=2)
    rand <- sample(seq(1, 6*length(diff),by=2), length(diff), replace=F)
    thresh <- .1
    negdiff <- diff*(diff< -thresh)
    posdiff <- diff*(diff>thresh)
    middiff <- diff*(diff <thresh & diff > -thresh)
    colors <- rgb(posdiff+.75, .75, -negdiff+.75, maxColorValue=2)
    text(diff, rand, model$vocab[words], cex=text.cex*scale, col=colors, family=family)
    
    #X-axis
    text(.75*min(diff), -length(diff), as.character(plabels[2]), 
         col=rgb(.5,.5,1.6,2,maxColorValue=2), cex=2, pos=1)
    text(.75*max(diff), -length(diff), as.character(plabels[1]), 
         col=rgb(1.6,.5,.5,2,maxColorValue=2), cex=2, pos=1)
    segments(min(diff),-.5*length(diff),max(diff),-.5*length(diff))
  }
  
  ###############
  # Theta Histogram plots
  ###############
  if(type=="hist") {
    #setup components of the plot size.
    N <- length(topics)
    root <- ceiling(sqrt(N))
    oldpar <- par(no.readonly=TRUE)
    par(mfrow=c(root,root), oma=c(0,.5,1.5,0), mar=c(2,2,4,1))
    
    if(labeltype!="custom"){
      lab <- apply(lab, 1, commas)
      if(!contentcov) lab <- lab[topics]
    }
    if(is.null(topic.names)) {
      topic.names <- sprintf("Topic %i:", topics)
    }
    lab <- sprintf("%s %s", topic.names, lab)
    
    for(i in 1:length(topics)){
      theta_median <- median(model$theta[,topics[i]])
      if(is.null(xlab)) xlab <- ""
      #Now call the histogram.  Note this is a bit kludgey but the default for xlim keeps us
      #  from being able to pass the NULL value.  Because it references something internally
      #  it would be difficult to make it one call.  So instead we just split on condition.
      if(is.null(xlim)) {
        hist(model$theta[,i], main=lab[i], xlab = xlab, ylab = "", ylim=ylim,...)
      } else {
        hist(model$theta[,i], main=lab[i], xlab = xlab, ylab = "", xlim=xlim, ylim=ylim,...)
      }
      abline(v = theta_median, col = "red", lwd = 1,lty=2)
    }
    
    if(is.null(main)) {
      title("Distribution of MAP Estimates of Document-Topic Proportions", outer=TRUE)
    } else {
      title(main, outer=TRUE)
    }
    par(oldpar)
  }
}

  
        
