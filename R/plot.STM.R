#' Functions for plotting STM objects
#' 
#' Produces one of four types of plots for an STM object.  The default option
#' \code{"summary"} prints topic words with their corpus frequency.
#' \code{"labels"} is for easy printing of tables of indicative words for each
#' topic.  \code{"perspectives"} depicts differences between two topics,
#' content covariates or combinations. \code{"hist"} creates a histogram of the
#' expected distribution of topic proportions across the documents.
#' 
#' The function can produce three types of plots which summarize an STM object
#' which is chosen by the argument \code{type}.  \code{summary} produces a plot
#' which displays the topics ordered by their expected frequency across the
#' corpus.  \code{labels} plots the top words selected according to the chosen
#' criteria for each selected topics.  \code{perspectives} plots two topic or
#' topic-covariate combinations.  Words are sized proportional to their use
#' within the plotted topic-covariate combinations and oriented along the
#' X-axis based on how much they favor one of the two configurations.  If the
#' words cluster on top of each other the user can either set the plot size to
#' be larger or shrink the total number of words on the plot.  The vertical
#' configuration of the words is random and thus can be rerun to produce
#' different results each time. Note that \code{perspectives} plots do 
#' not use any of the labeling options directly. \code{hist} plots a histogram of the MAP
#' estimates of the document-topic loadings across all documents.  The median
#' is also denoted by a dashed red line.
#' 
#' @param x Model output from stm.
#' @param type Sets the desired type of plot.  See details for more
#' information.
#' @param n Sets the number of words used to label each topic.  In perspective
#' plots it approximately sets the total number of words in the plot.  The
#' defaults are 3, 20 and 25 for \code{summary}, \code{labels} and
#' \code{perspectives} respectively.  n must be greater than or equal to 2
#' @param topics Vector of topics to display.  For plot perspectives this must
#' be a vector of length one or two. For the other two types it defaults to all
#' topics.
#' @param labeltype Determines which option of \code{"prob", "frex", "lift",
#' "score"} is used for choosing the most important words.  See
#' \code{\link{labelTopics}} for more detail.  Passing an argument to
#' \code{custom.labels} will override this. Note that this does not apply to
#' \code{perspectives} type which always uses highest probability words.
#' @param frexw If "frex" labeltype is used, this will be the frex weight.
#' @param main Title to the plot
#' @param xlim Range of the X-axis.
#' @param ylim Range of the Y-axis.
#' @param xlab Labels for the X-axis.  For perspective plots, use
#' \code{plabels} instead.
#' @param family The Font family.  Most of the time the user will not need to
#' specify this but if using other character sets can be useful see \link{par}.
#' @param width Sets the width in number of characters used for string wrapping
#' in type \code{"labels"}
#' @param covarlevels A vector of length one or length two which contains the
#' levels of the content covariate to be used in perspective plots.
#' @param plabels This option can be used to override the default labels in the
#' perspective plot that appear along the x-axis.  It should be a character
#' vector of length two which has the left hand side label first.
#' @param text.cex Controls the scaling constant on text size.
#' @param custom.labels A vector of custom labels if labeltype is equal to
#' "custom".
#' @param topic.names A vector of custom topic names.  Defaults to "Topic #: ".
#' @param \dots Additional parameters passed to plotting functions.
#' @seealso \code{\link{plotQuote}}, \code{\link{plot.topicCorr}}
#' @references Roberts, Margaret E., Brandon M. Stewart, Dustin Tingley,
#' Christopher Lucas, Jetson Leder-Luis, Shana Kushner Gadarian, Bethany
#' Albertson, and David G. Rand.  "Structural Topic Models for Open-Ended
#' Survey Responses." American Journal of Political Science 58, no 4 (2014):
#' 1064-1082.
#' @examples
#' \donttest{
#' 
#' #Examples with the Gadarian Data
#' plot(gadarianFit)
#' plot(gadarianFit,type="labels")
#' plot(gadarianFit, type="perspectives", topics=c(1,2))
#' plot(gadarianFit,type="hist")
#' }
#' @export plot.STM
#' @export
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
      text(frequency[invrank[i]]+ min(2*max(frequency), 1)/100, i , lab[invrank[i]],family=family,pos=4, cex=text.cex)
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
    colors <- grDevices::rgb(posdiff+.75, .75, -negdiff+.75, maxColorValue=2)
    text(diff, rand, model$vocab[words], cex=text.cex*scale, col=colors, family=family)
    
    #X-axis
    text(.75*min(diff), -length(diff), as.character(plabels[2]), 
         col=grDevices::rgb(.5,.5,1.6,2,maxColorValue=2), cex=2, pos=1)
    text(.75*max(diff), -length(diff), as.character(plabels[1]), 
         col=grDevices::rgb(1.6,.5,.5,2,maxColorValue=2), cex=2, pos=1)
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