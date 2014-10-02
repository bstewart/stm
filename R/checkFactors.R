  ################################################################################# 
# FACTOR CHECK
# Control for topics that load exclusively on a single (discrete) covariate level.
# Author: Antonio Coppola (acoppola@college.harvard.edu)
# June 17, 2014
################################################################################# 
# Inputs:
# - stmobject: The STM model object.
# - covariate: A vector containing the relevant document-level covariate. The 
#              length of this vector must be equal to the number of documents in
#              the model.
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

checkFactors <- function(stmobject, covariate=NULL, tolerance=0.01, 
                         reporting=0.99){

  # Check validity of tolerance argument.
  if (tolerance < 1e-6){
    stop("Tolerance value too low.")
  }

  # Initialize objects.
  K <- stmobject$settings$dim$K                       # Number of topics
  theta <- stmobject$theta                            # Theta matrix
  X <- covariate                                      # Covariate vector
  mylevels <- levels(factor(X))                         # Unique covariate levels
  indicators <- matrix(nrow=length(mylevels), ncol=K)   # Empty matrix
  proportions <- matrix(nrow=length(mylevels), ncol=K)  # Empty matrix
  err.tot <- 0                                        # Errors count

  # Control dimensions. The length of X must be equal to the # of documents.
  if (nrow(theta) != length(X))
    stop("Number of documents in model and covariate vector not corresponding.")

  # The theta matrix in the STM object (topical proportions by document) is 
  # subsetted by covariate levels. Indicators are computed specifying whether 
  # the entries in the subsetted matrix are higher than the tolerance value.
  indices <- apply(as.array(mylevels), 1, function(x) which(X == x))
  thetas.collapsed <- sapply(indices, function(x) theta[x,,drop=FALSE])

  for (k in seq(K)){
    for (i in seq(length(mylevels))){
      indicators[i, k] <- length(which(thetas.collapsed[[i]][ , k] > tolerance)) 
    }
  }

  # The indicators are aggregated by topic, summing over the columns of the 
  # subsetted matrix.
  tot.indicators <- colSums(indicators)

  # The aggregate values are divided by the total number of observations above 
  # the tolerance value in the corresponding column of the non-subsetted theta 
  # matrix to obtain the topical proportion by covariate level. If any of the
  # computed proportions are above the reporting threshold, an error is signaled.

  proportions <- sweep(indicators, MARGIN=2, STATS=tot.indicators, FUN="/") # Find errors by column
  errors <- t(apply(proportions, 1, function(x) x> reporting))
  err.tot <- length(which(errors))

  #If there are errors, print them
  if (err.tot == 0){
    cat("No abberrations found at the given reporting threshold.")
  }
  else{
    errorList <- which(errors, arr.ind = T)
    colnames(errorList) = c('Level', 'Topic')
    errorList[,1] <- mylevels[errorList[,1]]
    print('Topics/Levels Above Threshold:')
    print(errorList)
  }
  # Add row and column names.
  rownames(proportions) <- mylevels
  colnames(proportions) <- paste('Topic', seq(k))

  # Return output.
  output <- list(total.errors=err.tot, topical.proportions=proportions)
  return(invisible(output))
}