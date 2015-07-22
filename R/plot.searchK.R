plot.searchK<-function(x, ...){
  oldpar <- par(no.readonly=TRUE)
  g <- x$results
  par(mfrow=c(2,2),mar=c(4,4,4,4),oma=c(2,2,2,2))
  
  plot(g$K,g$heldout,type="p", main="Held-Out Likelihood", xlab="Number of Topics (K)", ylab="Held-Out Likelihood")
  lines(g$K,g$heldout,lty=1,col=1)

  plot(g$K,g$residual,type="p", main="Residuals", xlab="Number of Topics (K)", ylab="Residuals")
  lines(g$K,g$residual,lty=1,col=1 )
  
  
  plot(g$K,g$semcoh,type="p", main="Semantic Coherence", xlab="Number of Topics (K)", ylab="Semantic Coherence")
  lines(g$K,g$semcoh,lty=1,col=1 ) 
  
  #plot(g$K,g$exclus,type="n", main="Exclusivity", xlab="Number of Topics (K)", ylab="Exclusivity")
  #lines(g$K,g$exclus,lty=1,col=1 )  
  
  #plot(g$K,g$bound,type="n", main="Bound", xlab="Number of Topics (K)", ylab="Bound")
  #lines(g$K,g$bound,lty=1,col=1 ) 
  
  plot(g$K,g$lbound,type="p", main="Lower Bound", xlab="Number of Topics (K)", ylab="Lower Bound")
  lines(g$K,g$lbound,lty=1,col=1 ) 

  title("Diagnostic Values by Number of Topics", outer=TRUE)  
  par(oldpar)
}
