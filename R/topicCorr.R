topicCorr <- function(model, method=c("simple", "huge"), 
                      cutoff=.01, verbose = TRUE) {
  method <- match.arg(method)
  out <- list()
  if(method=="simple") {
    cormat <- cor(model$theta)
    adjmat <- ifelse(cormat > cutoff,1,0)
    out$posadj <- adjmat
    out$poscor <- cormat*adjmat
    out$cor <- ifelse(abs(cormat)> cutoff, cormat, 0)
  }
  if(method=="huge") {
    if(!requireNamespace("huge", quietly=TRUE)) stop("Install the huge package to use this function")
    X.npn <- huge::huge.npn(model$theta,verbose=verbose) # Nonparanormal
    out.npn <- huge::huge(X.npn,nlambda=30, verbose=verbose)
    ric.npn <- huge::huge.select(out.npn, verbose=verbose)
    MLE <- cor(model$theta)
    out$posadj <- ric.npn$refit*(MLE>0)
    out$poscor <- ric.npn$refit*(MLE>0)*MLE
    out$cor <- ric.npn$refit*MLE
  }
  class(out) <- "topicCorr"
  return(out)
}


plot.topicCorr <- function(x, topics=NULL,
                           vlabels=NULL, layout=NULL,
                           vertex.color="green", vertex.label.cex=.75, 
                           vertex.label.color="black", ...){
  if(!requireNamespace("igraph", quietly=TRUE)) stop("Install the igraph package to use this function.")
  if(is.null(topics)) topics <- 1:nrow(x$posadj)
  x <- x$posadj[topics, topics]
  
  g <- igraph::graph.adjacency(x, mode="undirected", weighted=TRUE, diag=FALSE)
  if(is.null(vlabels)) vlabels <-  paste("Topic", topics)
  igraph::E(g)$size <- 1
  igraph::E(g)$lty <- 2
  igraph::E(g)$color <- "black"
  igraph::V(g)$label <- vlabels
  if(is.null(layout)) layout <- igraph::layout.fruchterman.reingold
  igraph::plot.igraph(g, layout=layout, vertex.color=vertex.color, vertex.label.cex=vertex.label.cex, 
       vertex.label.color=vertex.label.color, ...)
}
