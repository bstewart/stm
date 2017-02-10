#' Summary Function for the STM objects
#' 
#' Function to report on the contents of STM objects
#' 
#' Summary prints a short statement about the model and then runs
#' \code{\link{labelTopics}}.
#' 
#' @aliases summary.STM print.STM
#' @param object An STM object.
#' @param \dots Additional arguments affecting the summary
#' @method summary STM
#' @export
summary.STM <- function(object,...) {
  toprint <- sprintf("A topic model with %i topics, %i documents and a %i word dictionary.\n", 
                     object$settings$dim$K, 
                     object$settings$dim$N, 
                     object$settings$dim$V)
  cat(toprint)
  labelTopics(object)
}

#' @method print STM
#' @export
print.STM <- function(x,...) {
  toprint <- sprintf("A topic model with %i topics, %i documents and a %i word dictionary.\n", 
                     x$settings$dim$K, 
                     x$settings$dim$N, 
                     x$settings$dim$V)
  cat(toprint)
}