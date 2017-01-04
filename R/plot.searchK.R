#' Plots diagnostic values resulting from searchK
#' 
#' Takes the result of searchK and produces a set of plots for evaluating
#' optimal topic numbers via visual representation of diagnostic functions.
#' 
#' 
#' @param x A searchK object, containing the diagnostic information of an stm
#' with a variety of topics.
#' @param ...  additional arguments for S3 compatability.
#' @examples
#' 
#' \dontrun{
#' 
#' K<-c(5,10,15) 
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' documents <- out$documents
#' vocab <- out$vocab
#' meta <- out$meta
#' set.seed(02138)
#' K<-c(5,10,15) 
#' kresult <- searchK(documents, vocab, K, prevalence=~treatment + s(pid_rep), data=meta)
#' 
#' plot(kresult)
#' }
#'  
#' @export
plot.searchK<-function(x, ...){
  oldpar <- par(no.readonly=TRUE)
  g <- x$results
  par(mfrow=c(2,2),mar=c(4,4,4,4),oma=c(2,2,2,2))
  
  plot(g$K,g$heldout,type="p", main="Held-Out Likelihood", xlab="Number of Topics (K)", ylab="Held-Out Likelihood")
  lines(g$K,g$heldout,lty=1,col=1)

  plot(g$K,g$residual,type="p", main="Residuals", xlab="Number of Topics (K)", ylab="Residuals")
  lines(g$K,g$residual,lty=1,col=1 )
  
  if(!is.null(g$semcoh)){
    plot(g$K,g$semcoh,type="p", main="Semantic Coherence", xlab="Number of Topics (K)", ylab="Semantic Coherence")
    lines(g$K,g$semcoh,lty=1,col=1 ) 
  }
  
  #plot(g$K,g$exclus,type="n", main="Exclusivity", xlab="Number of Topics (K)", ylab="Exclusivity")
  #lines(g$K,g$exclus,lty=1,col=1 )  
  
  #plot(g$K,g$bound,type="n", main="Bound", xlab="Number of Topics (K)", ylab="Bound")
  #lines(g$K,g$bound,lty=1,col=1 ) 
  
  plot(g$K,g$lbound,type="p", main="Lower Bound", xlab="Number of Topics (K)", ylab="Lower Bound")
  lines(g$K,g$lbound,lty=1,col=1 ) 

  title("Diagnostic Values by Number of Topics", outer=TRUE)  
  par(oldpar)
}