#' Plot some effects with loess
#' 
#' Plots a loess line of the topic proportions on a covariate inputted by the
#' user. This allows for a more flexible functional form for the relationship.
#' 
#' This function is considerably less developed than
#' \code{\link{plot.estimateEffect}} and we recommend using that function with
#' splines and high degrees of freedom where possible.  Computes standard
#' errors through the method of composition as in \code{\link{estimateEffect}}.
#' 
#' @param model An STM model object
#' @param topics Vector of topic numbers to plot by the covariate. E.g.,
#' c(1,2,3) would plot lines for topics 1,2,3.
#' @param covariate Covariate vector by which to plot topic proportions.
#' @param span loess span parameter.  See \code{\link{loess}}
#' @param level Desired coverage for confidence intervals
#' @param main Title of the plot, default is ""
#' @param xlab X-label, default is "Covariate"
#' @param ylab Y-label, default is "Topic Proportions"
#' @seealso \code{\link{plot.estimateEffect}}
#' @examples 
#' plotTopicLoess(gadarianFit, topics=1, covariate=gadarian$pid_rep)
#' @export
plotTopicLoess <- function(model, topics, covariate, span=1.5, level=.95,
                           main="", xlab="Covariate", ylab="Topic Proportions"){
  
  newdat <- data.frame(x=seq(min(covariate),max(covariate), length=100))
  draws <- vector(mode="list", length=length(topics))
  for(i in 1:100) {
    thetasims <- thetaPosterior(model, nsims=1, type="Global")
    thetasims <- do.call(rbind,thetasims)
    for(k in 1:length(topics)) {
      fit <- loess(thetasims[,k] ~ covariate, span=span)
      pred <- predict(fit, newdata=newdat$x, se=T)
      draws[[k]][[i]] <- sapply(1:length(pred$fit), function(x) rnorm(250,pred$fit[x],pred$se[x]))
    }
  }
  K <- length(topics)
  means <- matrix(0, nrow=K, ncol=100)
  low.ci <- matrix(0, nrow=K, ncol=100)
  upp.ci <- matrix(0, nrow=K, ncol=100)
  for(k in 1:length(topics)) {
    drawmat <- do.call(rbind, draws[[k]])
    means[k,] <- colMeans(drawmat)
    low.ci[k,] <- apply(drawmat, 2, quantile, (1-level)/2)
    upp.ci[k,] <- apply(drawmat, 2, quantile, 1- (1-level)/2)
  }

  plot(0, 0,col="white", ylab=ylab, xlab=xlab, 
       main=main, xlim=range(covariate), ylim=c(min(low.ci), max(upp.ci)))
  
  cols = grDevices::rainbow(length(topics))
  for(k in 1:K){
    lines(newdat$x, means[k,], col=cols[k])
    lines(newdat$x, low.ci[k,], col=cols[k], lty=2)
    lines(newdat$x, upp.ci[k,], col=cols[k], lty=2)
  }
}

.tjloesssim <- function(x, topic, covariate, newdat, span){
  fit <- loess(x[,topic] ~ covariate, span=span)
  pred <- predict(fit, newdata=newdat$x, se=T)
  return(list(fit=pred$fit, se=pred$se))
}

.TJciloess <- function(thetasims, topic, covariate, level=c(.025, .975), span=span, sims, newdat){
  betas <- lapply(thetasims, function (x) .tjloesssim(x, topic, covariate, newdat, span))
  betasims <- lapply(betas, function (x) apply(cbind(x[[1]],x[[2]]), 1, function (y) rnorm(sims, y[1], y[2])))
  cis <- sapply(1:sims, function (y) quantile(unlist(lapply(betasims, function (x) x[,y])), c(.025,.975)))
  mu <- sapply(1:sims, function (y) mean(unlist(lapply(betasims, function (x) x[,y]))))
  return(list(mean=mu, cis=cis))
}