plotModels <- function(models, xlab="Semantic Coherence", ylab="Exclusivity", labels=1:length(models$runout),...){
  if(length(models$runout[[1]]$beta$logbeta)<2){
    plot(0, 0, xlab=xlab, ylab=ylab, col="white", xlim=c(min(unlist(models$semcoh)), max(unlist(models$semcoh))),
         ylim=c(min(unlist(models$exclusivity)), max(unlist(models$exclusivity))),...)
    col = rainbow(length(models$runout))
    for(i in 1:length(models$runout)){
      points(models$semcoh[[i]], models$exclusivity[[i]], col=col[i], pch=16, cex=.75)
    }
    legend(min(unlist(models$semcoh)), max(unlist(models$exclusivity)),labels, col=col, pch=16)
    text(unlist(lapply(models$semcoh, mean)), unlist(lapply(models$exclusivity, mean)), labels, col=col)
  }
  if(length(models$runout[[1]]$beta$logbeta)>1){
    avgsc <- unlist(lapply(models$semcoh, mean))
    avgsp <- unlist(lapply(models$sparsity, mean))
    for(i in 1:length(models$runout)){
      print(paste("Model", i, "has on average", avgsc[i], "semantic coherence and", avgsp[i], "sparsity"))
    }
  }
}
