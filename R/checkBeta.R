checkBeta <- function(stmobject, tolerance=0.01){

  # Check validity of tolerance argument.
  if (tolerance < 1e-6){
    stop("Tolerance value too low.")
  }

  # Initialize objects
  betamatrix <- stmobject$beta$logbeta       # Log beta matrix (can be multiple)
  diagnostic <- list()                       # Empty list for working matrices
  topicErrorTotal <- list()                  # Count topics with errors
  problemTopics <- list()                    # Enumerate topics with errors
  wordErrorTotal <- list()                   # Count words with errors
  problemWords <- list()                     # Enumerate words with errors

  # Iterate through beta matrices and check for row abberrations
  for (i in seq(length(betamatrix))){
    listobject <- betamatrix[[i]]

    diagnostic[[i]] <- which(t(apply(listobject, 1, function(x) x > log(1 - tolerance))), arr.ind = TRUE)
    colnames(diagnostic[[i]]) <- c('Topic', 'Word')
    
    wordErrorTotal[[i]] <- nrow(diagnostic[[1]])
    problemWords[[i]]<-diagnostic[[i]][,c(2:1)]

    problemTopics[[i]] <- unlist(diagnostic[[i]][,1])
    topicErrorTotal[[i]] <- length(problemTopics[[i]])

    if(nrow(problemWords[[i]]) == 0){
      problemWords[[i]] <- rbind(problemWords[[i]], c('No Words', 'Load Exclusively Onto 1 Topic'))
    }
  }
  check <- any(!as.logical(wordErrorTotal))


  output <- list(problemTopics = problemTopics, topicErrorTotal = topicErrorTotal, problemWords = problemWords, wordErrorTotal = wordErrorTotal, checked = check)
  return(output)
}