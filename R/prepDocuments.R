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
  documents <- ijv.to.doc(triplet$i, triplet$j, triplet$v)
  docs.removed <- c()
  
  #Detect Missing Terms
  vocablist <- sort(unique(triplet$j))
  if(length(vocablist)>length(vocab)) {
    stop("Your documents object has more unique features than your vocabulary file has entries.")
  } 
  if(length(vocablist)<length(vocab)) {
    if(verbose) cat("Detected Missing Terms, renumbering \n")
    vocab <- vocab[vocablist]
    new.map <- cbind(vocablist, 1:length(vocablist))
    documents <- lapply(documents, function(d) {
      d[1,] <- new.map[match(d[1,], new.map[,1]),2]
      return(d)
    })
  }
  
  #Remove Words Appearing Only n Times
  wordcounts <- table(triplet$j)
  toremove <- which(wordcounts <= lower.thresh | wordcounts >= upper.thresh)
  keepers <- which(wordcounts > lower.thresh & wordcounts < upper.thresh)
  droppedwords <- vocab[toremove]
  if(length(toremove)) {
    if(verbose) cat("Removing Words Due to Frequency \n")
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
      if(verbose) cat("Removing Documents with No Words \n")
      documents <- documents[-docs.removed]
    }
    toprint <- sprintf("Your corpus now has %i documents and %i words.", length(documents), length(vocab))
    if(verbose) cat(toprint)
  }
  
  if(!is.null(docs.removed) & !is.null(meta)){
    meta<-meta[-docs.removed,]
  }
  #recast everything as an integer
  documents <- lapply(documents, function(x) matrix(as.integer(x), nrow=2))
  return(list(documents=documents, vocab=vocab, meta=meta, words.removed=droppedwords, docs.removed=docs.removed, wordcounts=wordcounts))
}