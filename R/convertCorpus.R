#' Convert \pkg{stm} formatted documents to another format
#' 
#' Takes an \pkg{stm} formatted documents and vocab object and returns formats
#' useable in other packages.
#' 
#' We also recommend the \pkg{quanteda} and \pkg{tm} packages for text preparation
#' etc.  The \code{convertCorpus} function is provided as a helpful utility for 
#' moving formats around, but if you intend to do text processing with a variety
#' of output formats, you likely want to start with \pkg{quanteda} or \pkg{tm}.
#'  
#' The various type conversions are described below: 
#' \describe{
#' \item{\code{type = "slam"}}{Converts to the simple triplet matrix
#' representation used by the \pkg{slam} package.  This is the format used
#' internally by \pkg{tm}.} 
#' \item{\code{type = "lda"}}{Converts to the format
#' used by the \pkg{lda} package.  This is a very minor change as the format in
#' \pkg{stm} is based on \pkg{lda}'s data representation.  The difference as
#' noted in \code{\link{stm}} involves how the numbers are indexed.
#' Accordingly this type returns a list containing the new documents object and
#' the unchanged vocab object.} 
#' \item{\code{type = "Matrix"}}{Converts to the
#' sparse matrix representation used by \pkg{Matrix}.  This is the format used
#' internally by numerous other text analysis packages.} }
#' 
#' If you want to write
#' out a file containing the sparse matrix representation popularized by David
#' Blei's \code{C} code \code{ldac} see the function \code{\link{writeLdac}}.
#' 
#' @param documents the documents object in \pkg{stm} format
#' @param vocab the vocab object in \pkg{stm} format
#' @param type the output type desired.  See Details.
#' @seealso \code{\link{writeLdac}} \code{\link{readCorpus}}
#' \code{\link{poliblog5k}}
#' @examples
#' #convert the poliblog5k data to slam package format
#' poliSlam <- convertCorpus(poliblog5k.docs, poliblog5k.voc, type="slam")
#' class(poliSlam)
#' poliMatrix <- convertCorpus(poliblog5k.docs, poliblog5k.voc, type="Matrix")
#' class(poliMatrix)
#' poliLDA <- convertCorpus(poliblog5k.docs, poliblog5k.voc, type="lda")
#' str(poliLDA)
#' @export
convertCorpus <- function(documents, vocab, type=c("slam", "lda", "Matrix")) {
  #a function to output our format to one of the existing frameworks
  type <- match.arg(type)
  switch(type,
         slam = convert.slam(documents,vocab),
         lda = convert.lda(documents,vocab),
         Matrix = convert.Matrix(documents, vocab))
}

convert.slam <- function(documents, vocab) {
  slam <- doc.to.ijv(documents)
  slam <- slam::simple_triplet_matrix(slam$i, slam$j, slam$v)
  colnames(slam) <- vocab
  return(slam)
}

convert.lda <- function(documents,vocab) {
  #just subtract one off from the indices
  documents<- lapply(documents, function(x) {
    x[1, ] <- as.integer(x[1, ] - 1)
    x})
  return(list(documents=documents, vocab=vocab))
}

convert.Matrix <- function(documents,vocab) {
  docs <- doc.to.ijv(documents)
  mat <- Matrix::sparseMatrix(docs$i,docs$j, x=docs$v)
  colnames(mat) <- vocab
  return(mat)
}