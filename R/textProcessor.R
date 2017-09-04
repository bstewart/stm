#' Process a vector of raw texts
#' 
#' Function that takes in a vector of raw texts (in a variety of languages) and
#' performs basic operations.  This function is essentially a wrapper \pkg{tm}
#' package where various user specified options can be selected.
#' 
#' This function is designed to provide a convenient and quick way to process a
#' relatively small volume texts for analysis with the package. It is designed
#' to quickly ingest data in a simple form like a spreadsheet where each
#' document sits in a single cell. If you have texts more complicated than a 
#' spreadsheet, we recommend you check out the excellent \pkg{readtext} package. 
#' 
#' The processor always strips extra white space but all other processing
#' options are optional.  Stemming uses the snowball stemmers and supports a
#' wide variety of languages.  Words in the vocabulary can be dropped due to
#' sparsity and stop word removal.  If a document no longer contains any words
#' it is dropped from the output.  Specifying meta-data is a convenient way to
#' make sure the appropriate rows are dropped from the corresponding metadata
#' file.
#' 
#' When the option \code{sparseLevel} is set to a number other than 1,
#' infrequently appearing words are removed.  When a term is removed from the
#' vocabulary a message will print to the screen (as long as \code{verbose} has
#' not been set to \code{FALSE}).  The message indicates the number of terms
#' removed (that is, the number of vocabulary entries) as well as the number of
#' tokens removed (appearences of individual words).  The function
#' \code{\link{prepDocuments}} provides additional methods to prune infrequent
#' words.  In general the functionality there should be preferred.
#' 
#' We emphasize that this function is a convenience wrapper around the
#' excellent \pkg{tm} package functionality without which it wouldn't be
#' possible.
#' 
#' @aliases textProcessor print.textProcessor head.textProcessor
#' summary.textProcessor
#' @param documents The documents to be processed.  A character vector where
#' each entry is the full text of a document (if passed as a different type
#' it will attempt to convert to a character vector).  
#' @param metadata Additional data about the documents.  Specifically a
#' \code{data.frame} or \code{matrix} object with number of rows equal to the
#' number of documents and one column per meta-data type. The column names are
#' used to label the metadata.  The metadata do not affect the text processing,
#' but providing the metadata object insures that if documents are dropped the
#' corresponding metadata rows are dropped as well.
#' @param lowercase Whether all words should be converted to lower case.
#' Defaults to TRUE.
#' @param removestopwords Whether stop words should be removed using the SMART
#' stopword list (in English) or the snowball stopword lists (for all other
#' languages). Defaults to TRUE.
#' @param removenumbers Whether numbers should be removed. Defaults to TRUE.
#' @param removepunctuation whether punctuation should be removed.  Defaults to
#' TRUE.
#' @param stem Whether or not to stem words. Defaults to TRUE
#' @param wordLengths From the \pkg{tm} package. An integer vector of length 2.
#' Words shorter than the minimum word length \code{wordLengths[1]} or longer
#' than the maximum word length \code{wordLengths[2]} are discarded. Defaults
#' to \code{c(3, Inf)}, i.e., a minimum word length of 3 characters.
#' @param sparselevel Removes terms where at least sparselevel proportion of
#' the entries are 0. Defaults to 1 which effectively turns the feature off.
#' @param language Language used for processing. Defaults to English. \code{tm}
#' uses the \code{SnowballC} stemmer which as of version 0.5 supports "danish
#' dutch english finnish french german hungarian italian norwegian portuguese
#' romanian russian spanish swedish turkish".  These can be specified as any on
#' of the above strings or by the three-letter ISO-639 codes.  You can also set
#' language to "na" if you want to leave it deliberately unspecified (see
#' documentation in \code{tm}) Note that languages listed here may not all have 
#' accompanying stopwords.  However if you have your own stopword list you can use
#' customstopwords below.
#' @param verbose If true prints information as it processes.
#' @param onlycharacter When TRUE, runs a regular expression substitution to
#' replace all non-alphanumeric characters. These characters can crash
#' textProcessor for some operating systems.  May remove foreign characters
#' depending on encoding. Defaults to FALSE.  Defaults to FALSE. Runs before
#' call to tm package.
#' @param striphtml When TRUE, runs a regular expression substitution to strip
#' html contained within <>.  Defaults to FALSE. Runs before call to tm
#' package.
#' @param customstopwords A character vector containing words to be removed.
#' Defaults to NULL which does not remove any additional words.  This function
#' is primarily for easy removal of application specific stopwords.  Note that
#' as with standard stopwords these are removed after converting everything to
#' lower case but before removing numbers, punctuation or stemming.  Thus words
#' to be removed should be all lower case but otherwise complete.
#' @param v1 A logical which defaults to \code{FALSE}.  If set to \code{TRUE} it
#' will use the ordering of operations we use used in Version 1.0 of the package.
#' @return \item{documents}{A list containing the documents in the stm format.}
#' \item{vocab }{Character vector of vocabulary.} \item{meta}{Data frame or
#' matrix containing the user-supplied metadata for the retained documents.}
#' @seealso \code{\link{readCorpus}}
#' @references Ingo Feinerer and Kurt Hornik (2013). tm: Text Mining Package. R
#' package version 0.5-9.1.
#' 
#' Ingo Feinerer, Kurt Hornik, and David Meyer (2008). Text Mining
#' Infrastructure in R. \emph{Journal of Statistical Software} 25(5): 1-54.
#' @examples
#' 
#' \dontrun{
#' 
#' head(gadarian)
#' #Process the data for analysis.
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' meta<-temp$meta
#' vocab<-temp$vocab
#' docs<-temp$documents
#' out <- prepDocuments(docs, vocab, meta)
#' docs<-out$documents
#' vocab<-out$vocab
#' meta <-out$meta
#' }
#' @export
textProcessor <- function(documents, metadata=NULL, 
                          lowercase=TRUE, removestopwords=TRUE, removenumbers=TRUE, removepunctuation=TRUE, stem=TRUE, 
                          wordLengths=c(3,Inf),sparselevel=1, language="en",
                          verbose=TRUE, onlycharacter=FALSE,striphtml=FALSE,
                          customstopwords=NULL, v1=FALSE) {
  if(!requireNamespace("tm",quietly=TRUE)) stop("Please install tm package to use this function. You will also need SnowballC if stemming.")
  if(!(utils::packageVersion("tm")>=0.6)) stop("Please install at least version 0.6 of the tm package.")
  if(stem) {
    if(!requireNamespace("SnowballC", quietly=TRUE)) stop("Please install SnowballC to use stemming.")
  }
  
  documents <- as.character(documents)
  
  if(striphtml){
    documents <- gsub('<.+?>', ' ', documents)
  }
  
  #remove non-visible characters
  documents <- stringr::str_replace_all(documents,"[^[:graph:]]", " ")
  
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
    if(utils::packageVersion("tm") >= "0.6") {
      txt <- tm::tm_map(txt, tm::content_transformer(tolower)) 
    } else {
      txt <- tm::tm_map(txt, tolower)
    }
  }
  
  if(!v1) {
    if(removepunctuation){
      if(verbose) cat("Removing punctuation... \n")
      txt <- tm::tm_map(txt, tm::removePunctuation, preserve_intra_word_dashes = TRUE) #Remove punctuation
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
  
  if(v1) {
    #return to the v1 style of removing punctuation right before stemming
    if(removepunctuation){
      if(verbose) cat("Removing punctuation... \n")
      txt <- tm::tm_map(txt, tm::removePunctuation, preserve_intra_word_dashes = TRUE) #Remove punctuation
    }    
  }
  
  if(stem){
    if(verbose) cat("Stemming... \n")
    txt <- tm::tm_map(txt, tm::stemDocument, language=language)
  }
  
  if(!is.null(metadata)) {
    for(i in 1:ncol(metadata)) {
       NLP::meta(txt, colnames(metadata)[i]) <- metadata[,i]
    }
  }
  
  #Make a matrix
  if(verbose) cat("Creating Output... \n")
  dtm <- tm::DocumentTermMatrix(txt, control=list(wordLengths=wordLengths))
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
    #if it is a type of data frame coerce it so we know its not a tibble or data.table
    if(inherits(metadata, "data.frame")) metadata <- as.data.frame(metadata)
    docindex <- unique(dtm$i)
    metadata <- NLP::meta(txt)[docindex, , drop = FALSE]
  }
  out <- read.slam(dtm) #using the read.slam() function in stm to convert
  
  ## It's possible that the processing has caused some documents to be
  ## dropped. These will be removed in the conversion from dtm to
  ## internal representation.  Better keep a record
  kept <- (1:length(documents) %in% unique(dtm$i))
  vocab <- as.character(out$vocab)
  out <- list(documents=out$documents, vocab=vocab, meta=metadata, docs.removed=which(!kept))
  class(out) <- "textProcessor"
  return(out)
}

#' @method print textProcessor
#' @export
print.textProcessor <- function(x,...) {
  toprint <- sprintf("A text corpus with %i documents, and an %i word dictionary.\n", 
                     length(x$documents), 
                     length(x$vocab))
  cat(toprint)
}

#' @method summary textProcessor
#' @export
summary.textProcessor <- function(object,...) {
  toprint <- sprintf("A text corpus with %i documents, and an %i word dictionary. Use str() to inspect object or see documentation \n", 
                     length(object$documents), 
                     length(object$vocab))
  cat(toprint)
}

#' @method head textProcessor
#' @export
head.textProcessor <- function(x,...) {
  for(i in 1:length(x)) {
    print(head(x[[i]]))
  }
}