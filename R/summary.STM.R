
summary.STM <- function(object,...) {
  toprint <- sprintf("A topic model with %i topics, %i documents and a %i word dictionary.\n", 
                     object$settings$dim$K, 
                     object$settings$dim$N, 
                     object$settings$dim$V)
  cat(toprint)
  labelTopics(object)
}

print.STM <- function(x,...) {
  toprint <- sprintf("A topic model with %i topics, %i documents and a %i word dictionary.\n", 
                     x$settings$dim$K, 
                     x$settings$dim$N, 
                     x$settings$dim$V)
  cat(toprint)
}