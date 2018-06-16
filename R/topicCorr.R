#' Estimate topic correlation
#' 
#' Estimates a graph of topic correlations using either a simple thresholding
#' measure or more sophisticated tests from the package \code{huge}.
#' 
#' We offer two estimation procedures for producing correlation graphs.  The
#' results of either method can be plotted using \code{\link{plot.topicCorr}}.
#' The first method is conceptually simpler and involves a simple thresholding
#' procedure on the estimated marginal topic proportion correlation matrix and
#' requires a human specified threshold.  The second method draws on recent
#' literature undirected graphical model estimation and is automatically tuned.
#' 
#' The \code{"simple"} method calculates the correlation of the MAP estimates
#' for the topic proportions \eqn{\theta} which yields the marginal correlation
#' of the mode of the variational distribution. Then we simply set to 0 those
#' edges where the correlation falls below the threshold.
#' 
#' An alternative strategy is to treat the problem as the recovery of edges in
#' a high-dimensional undirected graphical model. In these settings we assume
#' that observations come from a multivariate normal distribution with a sparse
#' precision matrix.  The goal is to infer which elements of the precision
#' matrix are non-zero corresponding to edges in a graph.  Meinhuasen and
#' Buhlmann (2006) showed that using sparse regression methods like the LASSO
#' it is possible to consistently identify edges even in very high dimensional
#' settings.
#' 
#' Selecting the option \code{"huge"} uses the \code{huge} package by Zhao and
#' Liu to estimate the graph.  We use a nonparanormal transformation of the
#' topic proportions (\eqn{\theta}) which uses semiparametric Gaussian copulas
#' to marginally transform the data.  This weakens the gaussian assumption of
#' the subsequent procedure.  We then estimate the graph using the Meinhuasen
#' and Buhlman procedure.  Model selection for the scale of the \eqn{L_1}
#' penalty is performed using the rotation information criterion (RIC) which
#' estimates the optimal degree of regularization by random rotations.  Zhao
#' and Lieu (2012) note that this selection approach has strong empirical
#' performance but is sensitive to under-selection of edges.  We choose this
#' metric as the default approach to model selection to reflect social
#' scientists' historically greater concern for false positive rates as opposed
#' to false negative rates.
#' 
#' We note that in models with low numbers of topics the simple procedure and
#' the more complex procedure will often yield identical results.  However, the
#' advantage of the more complex procedure is that it scales gracefully to
#' models with hundreds or even thousands of topics - specifically the set of
#' cases where some higher level structure like a correlation graph would be
#' the most useful.
#' 
#' @param model An STM object for which you want to estimate correlations
#' between topics.
#' @param method Method for estimating the graph.  \code{"simple"} simply
#' thresholds the covariances.  \code{"huge"} uses the semiparametric procedure
#' in the package \code{huge}.  See details below.
#' @param cutoff When using the simple method, this is the cutoff below which
#' correlations are truncated to zero.
#' @param verbose A logical which indicates whether information should be
#' printed to the screen when running \code{"huge"}.
#' @return \item{posadj}{K by K adjacency matrix where an edge represents
#' positive correlation selected by the model.} \item{poscor}{K by K
#' correlation matrix.  It takes values of zero where the correlation is either
#' negative or the edge is unselected by the model selection procedure.}
#' \item{cor}{K by K correlation matrix element-wise multiplied by the
#' adjacency matrix. Note that this will contain significant negative
#' correlations as well as positive correlations.}
#' @seealso \code{\link{plot.topicCorr}}
#' @references Lucas, Christopher, Richard A. Nielsen, Margaret E. Roberts,
#' Brandon M. Stewart, Alex Storer, and Dustin Tingley. "Computer-Assisted Text
#' Analysis for Comparative Politics." Political Analysis (2015).
#' 
#' T. Zhao and H. Liu. The huge Package for High-dimensional Undirected Graph
#' Estimation in R. Journal of Machine Learning Research, 2012
#' 
#' H. Liu, F. Han, M. Yuan, J. Lafferty and L. Wasserman. High Dimensional
#' Semiparametric Gaussian Copula Graphical Models. Annals of Statistics,2012
#' 
#' N. Meinshausen and P. Buhlmann. High-dimensional Graphs and Variable
#' Selection with the Lasso. The Annals of Statistics, 2006.
#' @export
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

#' Plot a topic correlation graph
#' 
#' Uses a topic correlation graph estimated by \code{\link{topicCorr}} and the
#' \code{igraph} package to plot a network where nodes are topics and edges
#' indicate a positive correlation.
#' 
#' Essentially a thin wrapper around the plotting functionality in the
#' \code{igraph} package. See package vignette for more details.
#' 
#' @param x A topicCorr model object.
#' @param topics A vector of topics to include in the plot, defaults to all.
#' @param vlabels A character vector of labels for the vertices.  Defaults to
#' "Topic #"
#' @param layout The layout algorithm passed to the \code{igraph} package.  It
#' will choose \code{layout.fruchterman.reingold} by default.  Note that to
#' pass an alternate algorithm you should load the \code{igraph} package first.
#' @param vertex.color Color of the vertices.
#' @param vertex.label.cex Controls the size of the labels.
#' @param vertex.label.color Controls the color of the labels.
#' @param vertex.size Controls the sizes of the vertices, either NULL, a scalar or a vector of the same length as number of topics.
#' @param \dots Additional parameters passed to \code{plot.graph.adjacency}
#' @seealso \code{\link{topicCorr}}
#' @references Csardi G, Nepusz T: The igraph software package for complex
#' network research, InterJournal, Complex Systems 1695. 2006.
#' http://igraph.sf.net
#' @examples
#' 
#' \dontrun{
#' 
#' #This function becomes more useful with larger numbers of topics.
#' #it is demonstrated here with a small model simply to show how the syntax works.
#' cormat <- topicCorr(gadarianFit)
#' plot(cormat)
#' }
#' @export plot.topicCorr
#' @export
plot.topicCorr <- function(x, topics=NULL,
                           vlabels=NULL, layout=NULL,
                           vertex.color="green", vertex.label.cex=.75, 
                           vertex.label.color="black",vertex.size=NULL, ...){
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
       vertex.label.color=vertex.label.color, vertex.size=vertex.size, ...)
}
