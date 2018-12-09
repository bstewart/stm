#' STM Corpus Coercion
#'
#' Convert a set of document term counts and associated metadata to
#' the form required for processing by the \code{\link{stm}} function.
#'
#' @param documents A documents-by-term matrix of counts, or a set of
#' counts in the format returned by \code{\link{prepDocuments}}. Supported
#' matrix formats include \pkg{quanteda} \link[quanteda]{dfm}
#' and \pkg{Matrix} sparse matrix objects in \code{"dgCMatrix"} or
#' \code{"dgTMatrix"} format.
#'
#' @param vocab Character vector specifying the words in the corpus in the
#' order of the vocab indices in documents. Each term in the vocabulary index
#' must appear at least once in the documents.  See \code{\link{prepDocuments}}
#' for dropping unused items in the vocabulary.  If \code{documents} is a
#' sparse matrix or \pkg{quanteda} \link[quanteda]{dfm} object, then \code{vocab} should not
#'  (and must not) be supplied.  It is contained already inside the column
#'  names of the matrix.
#'
#' @param data An optional data frame containing the prevalence and/or content
#' covariates.  If unspecified the variables are taken from the active
#' environment.
#'
#' @param \dots Additional arguments passed to or from other methods.
#'
#' @return A list with components \code{"documents"}, \code{"vocab"}, and
#' \code{"data"} in the form needed for further processing by the \code{stm}
#' function.
#'
#' @seealso \code{\link{prepDocuments}}, \code{\link{stm}}
#'
#' @examples
#' \dontrun{
#' library(quanteda)
#' gadarian_corpus <- corpus(gadarian, text_field = "open.ended.response")
#' gadarian_dfm <- dfm(gadarian_corpus, 
#'                      remove = stopwords("english"),
#'                      stem = TRUE)
#' asSTMCorpus(gadarian_dfm)
#' }
#' @export
asSTMCorpus <- function(documents, vocab, data = NULL, ...) {
  UseMethod("asSTMCorpus")
}

#' @method asSTMCorpus list
#' @export
#' @keywords internal
asSTMCorpus.list <- function(documents, vocab=NULL, data = NULL, ...) {
  list(documents = documents, vocab = vocab, data = data)
}

#' @method asSTMCorpus dfm
#' @export
#' @keywords internal
asSTMCorpus.dfm <- function(documents, vocab, data = NULL, ...) {
  if (!missing(vocab)) {
    # in case K was not specified by name, and it was confused with the
    # vocab argument (missing for dfm inputs)
    if (is.numeric(vocab) & length(vocab)==1) {
      stop("incorrect argument type for vocab, did you mean to specify K = ", vocab, "?")
    } else {
      stop("if documents is a dfm, do not specify vocab separately")
    }
  }

  # convert the dfm input as the first argument into the structure of the
  # older function where this is split into a list
  dfm_stm <- quanteda::convert(documents, to = "stm", docvars = data)
  if(is.null(data)) data <- dfm_stm[["meta"]]

  list(documents = dfm_stm[["documents"]], vocab = dfm_stm[["vocab"]],
       data = data)
}

#' @method asSTMCorpus dgCMatrix
#' @export
#' @importFrom methods as
#' @keywords internal
asSTMCorpus.dgCMatrix <- function(documents, vocab, data = NULL, ...) {
  if (!missing(vocab)) {
    # in case K was not specified by name, and it was confused with the
    # vocab argument (missing for dfm inputs)
    if (is.numeric(vocab) & length(vocab)==1) {
      stop("incorrect argument type for vocab, did you mean to specify K = ", vocab, "?")
    } else {
      stop("if documents is a matrix, do not specify vocab separately")
    }
  }

  # drop unused terms
  tot <- Matrix::colSums(documents)
  unseen <- which(tot == 0)
  if (length(unseen) > 0) {
      warning(sprintf("dropping %d unseen terms from the vocabulary", length(unseen)))
      documents <- documents[, -unseen, drop = FALSE]
  }

  # convert to docs-by-terms, one column per doc
  x <- as(Matrix::t(documents), "dgCMatrix")

  # get the vocab
  vocab <- rownames(x)

  # extract the term count information
  terms <- x@i + 1L
  counts <- as.integer(x@x)
  off <- x@p[-length(x@p)] + 1L
  len <- x@p[-1L] - x@p[-length(x@p)]

  # fill in the documents
  documents <- vector("list", ncol(x))
  names(documents) <- colnames(x)
  for (i in seq_along(documents)) {
      ix <- seq.int(off[[i]], length.out = len[[i]])
      documents[[i]] <- matrix(c(terms[ix], counts[ix]), nrow = 2L,
                               byrow = TRUE)
  }

  list(documents = documents, vocab = vocab, data = data)
}

#' @method asSTMCorpus dgTMatrix
#' @export
#' @keywords internal
asSTMCorpus.dgTMatrix <- function(documents, vocab, data = NULL, ...) {
    documents <- methods::as("dgCMatrix", documents)
    asSTMCorpus.dgCMatrix(documents, vocab, data, ...)
}
