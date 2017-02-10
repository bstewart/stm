#' Plots semantic coherence and exclusivity for high likelihood models
#' outputted from selectModel.
#' 
#' Plots semantic coherence and exclusivity for high likelihood models.  In the
#' case of models that include content covariates, prints semantic coherence
#' and sparsity.
#' 
#' Each model has semantic coherence and exclusivity values associated with
#' each topic.  In the default plot function, the small colored dots are
#' associated with a topic's semantic coherence and exclusivity.  Dots with the
#' same color as topics associated with the same model.  The average semantic
#' coherence and exclusivity is also plotted in the same color, but printed as
#' the model number associated with the output from selectModels().
#' 
#' With content covariates, the model does not output exclusivity because
#' exclusivity has been built in with the content covariates.  Instead, the
#' user should check to make sure that sparsity is high enough (typically
#' greater than .5), and then should select a model based on semantic
#' coherence.
#' 
#' @param models Output from selectModel.
#' @param labels Labels for each model.
#' @param xlab Character string that is x axis title. This will be semantic
#' coherence.
#' @param ylab Character string that is y axis title. This will be exclusivity.
#' @param ...  Other plotting parameters.
#' @export
plotModels <- function(models, xlab="Semantic Coherence", ylab="Exclusivity", labels=1:length(models$runout),...){
  if(!inherits(models, "selectModel")) {
    if(length(models)==1)   stop("plotModels only works for selectModel objects.")
    #we want to let it run if it was a case before selectModels had a class
    #thus just check if semantic coherence is in there.
    if(is.null(models$semcoh)) stop("plotModels only works for selectModel objects")
  }
  if(length(models$runout[[1]]$beta$logbeta)<2){
    plot(0, 0, xlab=xlab, ylab=ylab, col="white", xlim=c(min(unlist(models$semcoh)), max(unlist(models$semcoh))),
         ylim=c(min(unlist(models$exclusivity)), max(unlist(models$exclusivity))),...)
    col = grDevices::rainbow(length(models$runout))
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
