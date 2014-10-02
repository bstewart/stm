plot.searchK<-function(x, ...){
  g <- x$results
  par(mfrow=c(2,2),mar=c(4,4,4,4),oma=c(2,2,2,2))
  #jetson make better y labels.
  
  plot(g$K,g$heldout,type="n", main="Held-Out Likelihood", xlab="Number of Topics (K)", ylab="Held-Out Likelihood")
  lines(g$K,g$heldout,lty=1,col=1)

  plot(g$K,g$residual,type="n", main="Residuals", xlab="Number of Topics (K)", ylab="Residuals")
  lines(g$K,g$residual,lty=1,col=1 )
  
  
  #plot(g$K,g$semcoh,type="n", main="Semantic Coherence", xlab="Number of Topics (K)", ylab="Semantic Coherence")
  #lines(g$K,g$semcoh,lty=1,col=1 ) 
  
  #plot(g$K,g$exclus,type="n", main="Exclusivity", xlab="Number of Topics (K)", ylab="Exclusivity")
  #lines(g$K,g$exclus,lty=1,col=1 )  
  
  #plot(g$K,g$bound,type="n", main="Bound", xlab="Number of Topics (K)", ylab="Bound")
  #lines(g$K,g$bound,lty=1,col=1 ) 
  
  plot(g$K,g$lbound,type="n", main="L-Bound", xlab="Number of Topics (K)", ylab="L-Bound")
  lines(g$K,g$lbound,lty=1,col=1 ) 

  title("Diagnostic Values by Number of Topics", outer=TRUE)  
}
