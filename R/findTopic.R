#' Find topics that contain user specified words.
#' 
#' Find topics that contain user specified words.
#'
#' @param x The STM model object to be searched. May also be the output from
#' sageLabels.
#' @param list Character vector containing words to be searched.
#' @param n Number of words to consider
#' @param type Type of words to be searched.
#' @param verbose A logical indicating whether details should be printed to the
#' screen.
#' @seealso \code{\link{findThoughts}}
#' @examples
#' 
#' lab <- sageLabels(gadarianFit, n=5)
#' findTopic(lab, c("poor", "immigr", "peopl"))
#' findTopic(gadarianFit, c("poor", "immigr", "peopl"))
#' 
#' @export
findTopic <- function(x, list, n=20, type=c("prob", "frex", "lift","score"), verbose=TRUE) {
  type <- match.arg(type)
  if(inherits(x,"STM")) {
    x <- sageLabels(x, n=n)
  } else {
    if(!inherits(x,"sageLabels")) stop("x must be an STM or sageLabels object")
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