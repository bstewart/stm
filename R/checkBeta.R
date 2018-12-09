#' Looks for words that load exclusively onto a topic
#' 
#' Checks the log beta matrix for values too close to 0, which reflect words
#' that load onto a single topic.  
#' 
#' The function checks the log beta matrix for values that exceed the tolerance
#' threshold, indicating that a word has loaded onto a single topics. The output
#' gives the user lists of which topics have problems, which words in which
#' topics have problems, as well as a count of the total problems in topics and
#' the total number of problem words.
#' 
#' Note that if the tolerance value is below 1e-6, this function will throw an
#' error.
#' 
#' @param stmobject STM Model Output
#' @param tolerance User specified input reflecting closeness to 1.  E.g. a
#' tolerance of .01 will flag any values greater than .99.  Tolerance must be
#' above 1e-6.
#' @return \item{problemTopics}{A list of vectors, each vector corresponding to
#' the set of topics in the relevant beta matrix that contain words with too
#' high of a loading to that topic } \item{topicErrorTotal}{A list of integers,
#' each corresponding to the total number of topics with problems in the
#' relevant beta matrix} \item{problemWords}{A list of matrices, each
#' corresponding to a relevant beta matrix, which gives the topic and word
#' index of each word with too high of a topic loading} \item{wordErrorTotal}{A
#' list of integers, each corresponding to the total words with problems for
#' the relevant beta matrix} \item{check}{A boolean representing if the check
#' was passed. If wordErrorTotal is all 0s (no errors), check is True.}
#' @author Antonio Coppola
#' @examples 
#' checkBeta(gadarianFit)
#' @export
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