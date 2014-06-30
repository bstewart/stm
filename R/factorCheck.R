################################################################################# 
# FACTOR CHECK
# Control for topics that load exclusively onto a single covariate's factor.
# Author: Antonio Coppola (acoppola@college.harvard.edu)
# June 17, 2014
################################################################################# 
# Inputs:
# - stmobject: The STM model object.
# - covariate: A vector containing the relevant covariate levels.
# - tolerance: The minimum topical proportion (by document) to consider a document 
#              as containing a particular topic.
# - reporting: The minimum threshold for topical proportion (by covariate level)
#              to report an abberration.
################################################################################# 
# Outputs:
#   A list of two items:
#     - topical.proportions: A matrix reporting topical proportions by covariate
#                           levels. Rows correspond to levels, columns to topics.
#     - total.errors: The number of matrix entries above reporting threshold.
################################################################################# 

factorCheck <- function(stmobject, covariate=NULL, tolerance=0.01, 
                         reporting=0.99){
 
  # Check validity of tolerance argument
  if (tolerance < 1e-6){
    stop("Tolerance value too low.")
  }
  
  # Initialize objects
  K <- stmobject$settings$dim$K
  theta <- stmobject$theta
  X <- covariate
  levels <- levels(factor(X))
  indicators <- matrix(nrow=length(levels), ncol=K)
  proportions <- matrix(nrow=length(levels), ncol=K)
  err.tot <- 0
  
  # Control dimensions
  if (nrow(theta) != length(X))
    stop("Number of documents in model and covariate vector not corresponding.")
  
  # Loop through topics and covariate levels to find indicators  
  for (k in seq(K)){
    for (i in seq(length(levels))){
      indices <- which(X == levels[i])
      theta.collapsed <- theta[indices, ]
      indicators[i, k] <- length(which(theta.collapsed[ , k] > tolerance)) 
    }
  }
  
  # Find topical proportions by covariate level, and abberrations
  tot.indicators <- colSums(indicators)
  for (k in seq(K)){
    for (i in seq(length(levels))){
      proportions[i, k] <- indicators[i, k] / tot.indicators[k]
      if (proportions[i, k] > reporting){
        cat(paste("Topical proportion for covariate level ", 
                  levels[i], " and topic ", k, 
                  " is above reporting threshold.\n", sep=""))
        err.tot <- err.tot + 1
      }
    }
  }
  
  if (err.tot == 0)
    cat("No abberrations found at the given reporting threshold.")
  
  # Add row and column names
  rownames(proportions) <- levels
  cols <- c()
  length(cols) <- ncol(proportions)
  for (i in seq(ncol(proportions))){
    cols[i] <- paste("Topic", i)
  }
  colnames(proportions) <- cols
  
  output <- list(total.errors=err.tot, topical.proportions=proportions)
  return(invisible(output))
}