
#Takes a character vector with one entry per document and its metadata
textProcessor <- function(documents, metadata=NULL, 
                          lowercase=TRUE, removestopwords=TRUE, removenumbers=TRUE, removepunctuation=TRUE, stem=TRUE, 
                          sparselevel=1, language="en",
                          verbose=TRUE, onlycharacter=FALSE,striphtml=FALSE,
                          customstopwords=NULL) {
  if(!requireNamespace("tm",quietly=TRUE)) stop("Please install tm package to use this function. You will also need SnowballC if stemming.")
  if(stem) {
    if(!requireNamespace("SnowballC", quietly=TRUE)) stop("Please install SnowballC to use stemming.")
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
  
  if(striphtml){
  documents <- gsub('<.+?>', ' ', documents)
  }
  if(onlycharacter){
  documents <- gsub("[^[:alnum:]///' ]", " ", documents)
  }
  
  if(verbose) cat("Building corpus... \n")
  txt <- tm::VCorpus(tm::VectorSource(documents), readerControl=list(language= language))
  #Apply filters
  txt <- tm::tm_map(txt, tm::stripWhitespace)
  
  if(lowercase){
    if(verbose) cat("Converting to Lower Case... \n")
    #Convert to Lower case
    #(Note that this is slightly more complicated due to API change in tm)
    if(packageVersion("tm") >= "0.6") {
      txt <- tm::tm_map(txt, tm::content_transformer(tolower)) 
    } else {
      txt <- tm::tm_map(txt, tolower)
    }
  }
  if(removestopwords){
    if(verbose) cat("Removing stopwords... \n")
    txt <- tm::tm_map(txt, tm::removeWords, tm::stopwords(language)) #Remove stopwords
  }
  if(!is.null(customstopwords)) {
    if(verbose) cat("Remove Custom Stopwords...\n")
    txt <- tm::tm_map(txt, tm::removeWords, customstopwords)
  }
  if(removenumbers){
    if(verbose) cat("Removing numbers... \n")
    txt <- tm::tm_map(txt, tm::removeNumbers) #Remove numbers
  }
  if(removepunctuation){
    if(verbose) cat("Removing punctuation... \n")
    txt <- tm::tm_map(txt, tm::removePunctuation) #Remove punctuation
  }
  if(stem){
    if(verbose) cat("Stemming... \n")
    txt <- tm::tm_map(txt, tm::stemDocument, language=language)
  }
  
  if(!is.null(metadata)) {
    for(i in 1:ncol(metadata)) {
     if(packageVersion("tm") >= "0.6") {   
       NLP::meta(txt, colnames(metadata)[i]) <- metadata[,i]
     } else {
       tm::meta(txt, colnames(metadata)[i]) <- metadata[,i]
     }
    }
  }
  
  #Make a matrix
  if(verbose) cat("Creating Output... \n")
  dtm <- tm::DocumentTermMatrix(txt)
  if(sparselevel!=1) {
    ntokens <- sum(dtm$v)
    V <- ncol(dtm)
    dtm <- tm::removeSparseTerms(dtm, sparselevel) #remove terms that are sparse
    if(ncol(dtm) < V & verbose) {
      message <- sprintf("Removed %i of %i terms (%i of %i tokens) due to sparselevel of %s \n", 
                       V-ncol(dtm), V,
                       ntokens-sum(dtm$v), ntokens,
                       sparselevel)
      cat(message)
    }
  }
  #If there is metadata we need to remove some documents
  if(!is.null(metadata)) {
    docindex <- unique(dtm$i)
    if(packageVersion("tm") >= "0.6") {
      metadata <- NLP::meta(txt)[docindex,]
    } else {
      metadata <- tm::meta(txt)[docindex,]
    }
  }
  out <- read.slam(dtm) #using the read.slam() function in stm to convert
  vocab <- as.character(out$vocab)
  return(list(documents=out$documents, vocab=vocab, meta=metadata))
}
