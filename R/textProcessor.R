
#Takes a character vector with one entry per document and its metadata
textProcessor <- function(documents, metadata=NULL, 
                          lowercase=TRUE, removestopwords=TRUE, removenumbers=TRUE, removepunctuation=TRUE, stem=TRUE, 
                          sparselevel=1, language="en",
                          verbose=TRUE) {
  if(!require(tm,quietly=TRUE)) stop("Please install tm package to use this function. You will also need SnowballC if stemming.")
  if(stem) {
    if(!require(SnowballC, quietly=TRUE)) stop("Please install SnowballC to use stemming.")
  }
  
  #If there is only one item assume its a url and load it.
  if(length(documents)==1) {
    filelist <- list.files(path=documents, full.names=TRUE, recursive=TRUE)
    documents <- vector(length=length(filelist))
    if(verbose) cat(sprintf("Loading %i files from directory...\n", length(documents)))
    for(i in 1:length(filelist)) {
      documents[i] <- paste(readLines(filelist[i]), collapse=" ")
    } 
  } else {
    documents <- as.character(documents)
  }
  
  if(verbose) cat("Building corpus... \n")
  txt <- VCorpus(VectorSource(documents), readerControl=list(language= language))
  #Apply filters
  txt <- tm_map(txt, stripWhitespace)
  
  if(lowercase){
    if(verbose) cat("Converting to Lower Case... \n")
    #Convert to Lower case
    #(Note that this is slightly more complicated due to API change in tm)
    if(packageVersion("tm") >= "0.6") {
      txt <- tm_map(txt, content_transformer(tolower)) 
    } else {
      txt <- tm_map(txt, tolower)
    }
  }
  if(removestopwords){
    if(verbose) cat("Removing stopwords... \n")
    txt <- tm_map(txt, removeWords, stopwords(language)) #Remove stopwords
  }
  if(removenumbers){
    if(verbose) cat("Removing numbers... \n")
    txt <- tm_map(txt, removeNumbers) #Remove numbers
  }
  if(removepunctuation){
    if(verbose) cat("Removing punctuation... \n")
    txt <- tm_map(txt, removePunctuation) #Remove punctuation
  }
  if(stem){
    if(verbose) cat("Stemming... \n")
    txt <- tm_map(txt, stemDocument, language=language)
  }
  
  if(!is.null(metadata)) {
    for(i in 1:ncol(metadata)) {
      meta(txt, colnames(metadata)[i]) <- metadata[,i]
    }
  }
  
  #Make a matrix
  if(verbose) cat("Creating Output... \n")
  dtm <- DocumentTermMatrix(txt)
  if(sparselevel!=1) {
    V <- ncol(dtm)
    dtm <- removeSparseTerms(dtm, sparselevel) #remove terms that are sparse
    if(ncol(dtm) < V & verbose) {
      message <- sprintf("Removed %i of %i words due to sparselevel of %s \n", 
                       V-ncol(dtm), V, sparselevel)
      cat(message)
    }
  }
  #If there is metadata we need to remove some documents
  if(!is.null(metadata)) {
    docindex <- unique(dtm$i)
    metadata <- meta(txt)[docindex,]
  }
  out <- read.slam(dtm) #using the read.slam() function in stm to convert
  vocab <- as.character(out$vocab)
  return(list(documents=out$documents, vocab=vocab, meta=metadata))
}
