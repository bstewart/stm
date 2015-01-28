findTopic <- function(x, list, n=20, type=c("prob", "frex", "lift","score"), verbose=TRUE) {
  type <- match.arg(type)
  if(class(x)=="STM") {
    x <- sageLabels(x, n=n)
  } else {
    if(class(x)!="sageLabels") stop("x must be an STM or sageLabels object")
  }
  counts <- apply(x$marginal[[type]],1, function(w) sum(list%in%w))
  if(max(counts)==0) {
    if(verbose) cat("No topics contained any words in the list.")
    return(invisible(0))
  } else {
    index <- which(counts==max(counts))
    if(verbose) cat(sprintf("%i topics contained %i words in the list: %s", length(index),
                max(counts), paste(index, collapse=", ")))
    return(invisible(index))
  }
}
