#' Write a .ldac file
#' 
#' A function for writing documents out to a .ldac formatted file.
#' 
#' This is a simple convenience function for writing out document corpora.
#' Files can be read back into R using \code{\link{readCorpus}} or simply used
#' for other programs.  The output is a file in the \code{.ldac} sparse matrix
#' format popularized by Dave Blei's C code for LDA.
#' 
#' @param documents A documents object to be written out to a file.  Object
#' must be a list of with each element corresponding to a document.  Each
#' document is represented as an integer matrix with two rows, and columns
#' equal to the number of unique vocabulary words in the document.  The first
#' row contains the 1-indexed vocabulary entry and the second row contains the
#' number of times that term appears
#' @param file A character string giving the name of the file to be written.
#' This object is passed directly to the argument \code{con} in
#' \code{\link{writeLines}} and thus can be a connection object as well.
#' @param zeroindex If \code{TRUE} (the default) it subtracts one
#' from each index.  If \code{FALSE} it uses the indices as given.  The
#' standard \code{.ldac} format indexes from 0 as per standard convention in
#' most languages.  Our documents format indexes from 1 to abide by conventions
#' in \code{R}.  This option converts to the zero index by default.
#' @seealso \code{\link{readCorpus}}
#' @examples
#' 
#' \dontrun{
#' 
#' #Convert the gadarian data into documents format
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' documents <- temp$documents
#' #Now write out to an ldac file
#' writeLdac(documents, file="gadarian.ldac")
#' }
#' @export
writeLdac <- function(documents, file, zeroindex=TRUE) {
  if(missing(file)) stop("Please specify a filename with argument file.")
  lines <- lapply(documents, function(x) {
    Nd <-ncol(x)
    if(zeroindex) {
      paste(Nd, paste(x[1,]-1, x[2,], sep=":", collapse=" "))
    } else {
      paste(Nd, paste(x[1,], x[2,], sep=":", collapse=" "))
    }
  })
  lines <- unlist(lines)
  writeLines(lines, con=file)
}