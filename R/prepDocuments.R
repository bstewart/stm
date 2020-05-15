#' Prepare documents for analysis with \code{stm}
#' 
#' Performs several corpus manipulations including removing words and
#' renumbering word indices (to correct for zero-indexing and/or unused words
#' in the vocab vector).
#' 
#' The default setting \code{lower.thresh=1} means that words which appear in
#' only one document will be dropped.  This is often advantageous as there is
#' little information about these words but the added cost of including them in
#' the model can be quite large.  In many cases it will be helpful to set this
#' threshold considerably higher.  If the vocabulary is in excess of 5000
#' entries inference can slow quite a bit.
#' 
#' If words are removed, the function returns a vector of the original indices
#' for the dropped items.  If it removed documents it returns a vector of doc
#' indices removed. Users with accompanying metadata or texts may want to drop
#' those rows from the corresponding objects.
#' 
#' The behavior is such that when \code{prepDocuments} drops documents their
#' corresponding rows are deleted and the row names are not renumbered.  We however
#' do not recommend using rownames for joins- instead the best practice is to either
#' keep a unique identifier in the \code{meta} object for doing joins or use something
#' like \pkg{quanteda} which has a more robust interface for manipulating the corpus
#' itself.
#' 
#' If you have any documents which are of length 0 in your original object the
#' function will throw an error. These should be removed before running the
#' function although please be sure to remove the corresponding rows in the
#' meta data file if you have one.  You can quickly identify the documents
#' using the code: \code{which(unlist(lapply(documents, length))==0)}.
#' 
#' @param documents List of documents. For more on the format see
#' \code{\link{stm}}.
#' @param vocab Character vector of words in the vocabulary.
#' @param meta Document metadata.
#' @param lower.thresh Words which do not appear in a number of documents
#' greater than lower.thresh will be dropped and both the documents and vocab
#' files will be renumbered accordingly.  If this causes all words within a
#' document to be dropped, a message will print to the screen at it will also
#' return vector of the documents removed so you can update your meta data as
#' well. See details below.
#' @param upper.thresh As with lower.thresh but this provides an upper bound.
#' Words which appear in at least this number of documents will be dropped.
#' Defaults to \code{Inf} which does no filtering.
#' @param subsample If an integer will randomly subsample (without replacement)
#' the given number of documents from the total corpus before any processing.
#' Defaults to \code{NULL} which provides no subsampling.  Note that the output
#' may have fewer than the number of requested documents if additional
#' processing causes some of those documents to be dropped.
#' @param verbose A logical indicating whether or not to print details to the
#' screen.
#' @return A list containing a new documents and vocab object.
#' \item{documents}{The new documents object for use with \code{stm}}
#' \item{vocab}{The new vocab object for use with \code{stm}} \item{meta}{The
#' new meta data object for use with \code{stm}. Will be the same if no
#' documents are removed.} \item{words.removed}{A set of indices corresponding
#' to the positions in the original vocab object of words which have been
#' removed.} \item{docs.removed}{A set of indices corresponding to the
#' positions in the original documents object of documents which no longer
#' contained any words after dropping terms from the vocab.}
#' \item{tokens.removed}{An integer corresponding to the number of unique
#' tokens removed from the corpus.} \item{wordcounts}{A table giving the the
#' number of documents that each word is found in of the original document set,
#' prior to any removal. This can be passed through a histogram for visual
#' inspection.}
#' @seealso \code{\link{plotRemoved}}
#' @examples
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' meta<-temp$meta
#' vocab<-temp$vocab
#' docs<-temp$documents
#' out <- prepDocuments(docs, vocab, meta)
#' @export
prepDocuments <- function(documents, vocab, meta=NULL, 
                           lower.thresh=1, upper.thresh=Inf, 
                           subsample=NULL,
                           verbose=TRUE) {
  #Functions:
  # 1) Detect and renumber zero-indexed data.
  # 2) Detect and renumber missing terms
  # 3) Remove words appearing only in [lower.thresh] documents.
  # 4) Remove words appear in upper.thresh or more of the documents.
  
  # Can also optionally subsample the data.
  
  #error check for inputs
  if((is.null(documents)|is.null(vocab))){
    stop("One of your file inputs has no data")
  }
  if(!inherits(documents,"list")) {
    stop("documents must be a list in stm() format.  See ?stm() for format.  
          See ?readCorpus() for tools for converting from popular formats")
  }
  
  if(!is.null(subsample)) {
    index <- sample(1:length(documents), subsample)
    documents <- documents[index]
    if(!is.null(meta)) meta <- meta[index, , drop = FALSE] 
  }
  
  #check that there are no 0 length documents
  len <- unlist(lapply(documents, length))
  if(any(len==0)) {
    stop("Some documents have 0 length.  Please check input. 
          See ?prepDocuments() for more info.")
  } 
   
  triplet <- doc.to.ijv(documents) #this also fixes the zero indexing.
  nms <- names(documents) 
  documents <- ijv.to.doc(triplet$i, triplet$j, triplet$v)
  names(documents) <- nms
  docs.removed <- c()
  
  #Detect Missing Terms
  miss.vocab <- NULL
  vocablist <- sort(unique(triplet$j))
  wordcounts <- tabulate(triplet$j)
  if(length(vocablist)>length(vocab)) {
    stop("Your documents object has more unique features than your vocabulary file has entries.")
  } 
  if(length(vocablist)<length(vocab)) {
    if(verbose) cat("Detected Missing Terms, renumbering \n")
    miss.vocab <- vocab[-vocablist]
    vocab <- vocab[vocablist]
    new.map <- cbind(vocablist, 1:length(vocablist))
    documents <- lapply(documents, function(d) {
      nm <- names(d)
      d[1,] <- new.map[match(d[1,], new.map[,1]),2]
      names(d) <- nm
      return(d)
    })
    wordcounts <- wordcounts[vocablist]
  }
  
  #Remove Words Appearing Only n Times
  toremove <- which(wordcounts <= lower.thresh | wordcounts >= upper.thresh)
  keepers <- which(wordcounts > lower.thresh & wordcounts < upper.thresh)
  droppedwords <- c(miss.vocab,vocab[toremove])
  if(length(toremove)) {
    if(verbose) cat(sprintf("Removing %i of %i terms (%i of %i tokens) due to frequency \n", 
                            length(toremove), length(wordcounts), sum(wordcounts[toremove]), sum(wordcounts)))
    vocab <- vocab[-toremove]
    remap <- 1:length(keepers)
    for(i in 1:length(documents)) {
      doc <- documents[[i]]
      dockeep <- doc[1,]%in%keepers
      doc <- doc[,dockeep,drop=FALSE]
      doc[1,] <- remap[match(doc[1,], keepers)]
      documents[[i]] <- doc
      if(ncol(doc)==0) docs.removed <- c(docs.removed,i)
    }
    if(length(docs.removed)) {
      if(verbose) cat(sprintf("Removing %i Documents with No Words \n", length(docs.removed)))
      documents <- documents[-docs.removed]
    }
    toprint <- sprintf("Your corpus now has %i documents, %i terms and %i tokens.", 
                       length(documents), length(vocab), sum(wordcounts[keepers]))
    if(verbose) cat(toprint)
  }
  
  if(!is.null(docs.removed) & !is.null(meta)){
    meta<-meta[-docs.removed, , drop = FALSE]
  }
  #recast everything as an integer
  documents <- lapply(documents, function(x) matrix(as.integer(x), nrow=2))
  return(list(documents=documents, vocab=vocab, meta=meta, 
              words.removed=droppedwords, docs.removed=docs.removed, 
              tokens.removed=sum(wordcounts[toremove]), wordcounts=wordcounts))
}