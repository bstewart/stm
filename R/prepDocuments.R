###
# A function for pre-processing documents
###
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
  if(class(documents)!="list") {
    stop("documents must be a list in stm() format.  See ?stm() for format.  
          See ?readCorpus() for tools for converting from popular formats")
  }
  
  if(!is.null(subsample)) {
    index <- sample(1:length(documents), subsample)
    documents <- documents[index]
    if(!is.null(meta)) meta <- meta[index,] 
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
    meta<-meta[-docs.removed,]
  }
  #recast everything as an integer
  documents <- lapply(documents, function(x) matrix(as.integer(x), nrow=2))
  return(list(documents=documents, vocab=vocab, meta=meta, 
              words.removed=droppedwords, docs.removed=docs.removed, 
              tokens.removed=sum(wordcounts[toremove]), wordcounts=wordcounts))
}