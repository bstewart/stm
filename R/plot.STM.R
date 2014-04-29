#Functions for plotting STM objects

plot.STM <- function(x, 
                     type=c("summary", "labels", "perspectives"),
                     n=NULL, topics=NULL,
                     labeltype="prob", frexw=.5,
                     main=NULL, xlim=NULL, ylim=NULL, xlab=NULL, family="",
                     width=80, 
                     covarlevels=NULL, plabels=NULL, text.cex=1, custom.labels=NULL,
                     ...){
  model <- x
  type <- match.arg(type)
  if(is.null(n)) n <- switch(type, 
                             summary=3, 
                             labels=20,
                             perspectives=25)
  contentcov <- length(model$beta$logbeta)!=1
  
  if(labeltype!="custom"){
    if(type %in% c("summary", "labels")) {
      if(is.null(topics)) topics <- 1:model$settings$dim$K
                                        #compute the labels
      lab <- labelTopics(model, topics=topics, n = n, frexweight=frexw)
      if(contentcov) {
        lab <- lab$topics
      } else {
        lab <- lab[[labeltype]]
      }    
    }
  }
  
  if(labeltype=="custom"){
    lab <- custom.labels
  }
  
  # Summary Plot
  # note at one point there was a threshold metric for choosing the number of words but
  # it was specific to highest frequency words so I opted to just remove it.
  if(type=="summary") {
    if(labeltype!="custom"){
      lab <- apply(lab, 1, commas)
      lab <- sprintf("Topic %i: %s", 1:length(lab), lab)
      lab <- lab[topics]
    }
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
      text(frequency[invrank[i]]+.05, i , lab[invrank[i]],family=family,pos=4, cex=text.cex)
    }
  }
  
  # Labels Plot
  if(type=="labels") {
    if(contentcov) stop("labels plot not yet implemented with content covariates.  See labelTopics for labels.")
    plot(c(0,0), type="n", main=main,
         ylim=c(1,length(topics)+.9), 
         xlim=c(1,n+1), 
         xaxt="n", yaxt="n", xlab="", ylab="",...)
    if(labeltype!="custom") lab <- apply(lab, 1, commas)
    for(i in 1:length(topics)){
      if(i!=length(topics)) lines(c(-1,n+3), c(i+1,i+1), lty=2)
      string <- paste0(strwrap(lab[rev(topics)[i]], width=width),collapse="\n")
      text(n/2, i + .5, 
           paste("Topic", rev(topics)[i], ":  \n", string), family=family, cex=text.cex)
    }
  }
  
  # Perspectives 
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
    plot(c(0,0), xlim=c(min(diff) - .1, max(diff)+.1), ylim=c(-2*length(diff),6*length(diff)), type="n", 
         xaxt="n", xlab="", main=main, yaxt="n", ylab="", bty="n",...)
    segments(0,0, 0, -2*length(diff),lty=2)
    rand <- sample(seq(1, 6*length(diff),by=2), length(diff), replace=F)
    thresh <- .1
    negdiff <- diff*(diff< -thresh)
    posdiff <- diff*(diff>thresh)
    middiff <- diff*(diff <thresh & diff > -thresh)
    colors <- rgb(posdiff+.75, .75, -negdiff+.75, maxColorValue=2)
    text(diff, rand, model$vocab[words], cex=scale, col=colors, family=family)
    
    #X-axis
    text(.75*min(diff), -length(diff), as.character(plabels[2]), 
         col=rgb(.5,.5,1.6,2,maxColorValue=2), cex=2, pos=1)
    text(.75*max(diff), -length(diff), as.character(plabels[1]), 
         col=rgb(1.6,.5,.5,2,maxColorValue=2), cex=2, pos=1)
    segments(min(diff),-.5*length(diff),max(diff),-.5*length(diff))
  }
}

  
        
