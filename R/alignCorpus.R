#' Align the vocabulary of a new corpus to an old corpus
#' 
#' Function that takes in a list of documents, vocab and (optionally) metadata
#' for a corpus of previously unseen documents and aligns them to an old vocabulary. 
#' Helps preprocess documents for \code{\link{fitNewDocuments}}.
#' 
#' When estimating topic proportions for previously unseen documents using
#' \code{\link{fitNewDocuments}} the new documents must have the same vocabulary
#' ordered in the same was as the original model.  This function helps with that
#' process.
#' 
#' Note: the code is not really built for speed or memory efficiency- if you are trying
#' to do this with a really large corpus of new texts you might consider building the object
#' yourself using \pkg{quanteda} or some other option.
#' 
#' @param new a list (such as those produced by \code{textProcessor} or \code{prepDocuments}) 
#' containing a list of documents in \code{\link{stm}} format, a character vector 
#' containing the vocabulary and optional a \code{data.frame} containing meta data.
#' These should be labelled \code{documents}, \code{vocab},and \code{meta} respectively.
#' This is the new set of unseen documents which will be returned with the vocab renumbered
#' and all words not appearing in \code{old} removed.
#' @param old.vocab a character vector containing the vocabulary that you want to align to.
#' In general this will be the vocab used in your original stm model fit which from an stm
#' object called \code{mod} can be accessed as \code{mod$vocab}.  
#' @param verbose a logical indicating whether information about the new corpus should be
#' printed to the screen. Defaults to \code{TRUE}.  
#' @return \item{documents}{A list containing the documents in the stm format.}
#' \item{vocab }{Character vector of vocabulary.} \item{meta}{Data frame or
#' matrix containing the user-supplied metadata for the retained documents.}
#' \item{docs.removed}{document indices (corresponding to the original data passed) of
#' documents removed because they contain no words}
#' \item{words.removed}{words dropped from \code{new}}
#' \item{tokens.removed}{the total number of tokens dropped from the new documents.}
#' \item{wordcounts}{counts of times the old vocab appears in the new documents}
#' \item{prop.overlap}{length two vector used to populate the message printed by verbose.}
#' @seealso \code{\link{prepDocuments}} \code{\link{fitNewDocuments}}
#' @examples
#' #we process an original set that is just the first 100 documents
#' temp<-textProcessor(documents=gadarian$open.ended.response[1:100],metadata=gadarian[1:100,])
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' set.seed(02138)
#' #Maximum EM its is set low to make this run fast, run models to convergence!
#' mod.out <- stm(out$documents, out$vocab, 3, prevalence=~treatment + s(pid_rep), 
#'               data=out$meta, max.em.its=5)
#' #now we process the remaining documents
#' temp<-textProcessor(documents=gadarian$open.ended.response[101:nrow(gadarian)],
#'                     metadata=gadarian[101:nrow(gadarian),])
#' #note we don't run prepCorpus here because we don't want to drop any words- we want 
#' #every word that showed up in the old documents.
#' newdocs <- alignCorpus(new=temp, old.vocab=mod.out$vocab)
#' #we get some helpful feedback on what has been retained and lost in the print out.
#' #and now we can fit our new held-out documents
#' fitNewDocuments(model=mod.out, documents=newdocs$documents, newData=newdocs$meta,
#'                 origData=out$meta, prevalence=~treatment + s(pid_rep),
#'                 prevalencePrior="Covariate")
#' @export
alignCorpus <- function(new, old.vocab, verbose=TRUE) {
  if(!is.character(old.vocab)) stop("old.vocab must be a character vector containing the original vocabulary.")
  if(!is.list(new)) stop("new must be a list containing at least a documents and vocab object.")
  #a function to take a new corpus and map it into the space defined by old
  id <- match(new$vocab, old.vocab, nomatch=0)
  num_notmatched <- sum(id==0)
  new_vocab_size <- length(new$vocab)
  
  #total tokens
  total_tokens <- sum(unlist(lapply(new$documents, function(x) sum(x[2,]))))
  
  #create storage for word counts
  wcts <- rep(0, length(old.vocab))
  #drop will store documents that lose all of their words
  drop <- rep(0, length(new$documents))
  for(i in 1:length(new$documents)) {
    doc <- new$documents[[i]]
    doc[1,] <- id[doc[1,]] #replace with new index
    #drop vocab that failed to match
    doc <- doc[,which(doc[1,] != 0), drop=FALSE]
    
    if(ncol(doc)==0) {
      drop[i] <- 1L
    } else {
      wcts[doc[1,]] <- wcts[doc[1,]] + doc[2,]
    }
    new$documents[[i]] <- doc
  }
  
  #have to figure out words removed before copying new vocab
  new$words.removed <- new$vocab[id==0]
  
  #drop documents that don't have any words
  if(any(drop==1)) {
    new$documents <- new$documents[drop==0]
    new$meta <- new$meta[drop==0,]
    new$vocab <- old.vocab
  }
  
  new_tokens <- sum(unlist(lapply(new$documents, function(x) sum(x[2,]))))
  
  #write out a bunch of nice things
  new$docs.removed <- which(drop==1)
  new$tokens.removed  <- total_tokens - new_tokens
  new$wordcounts <- wcts 
  new$prop.overlap <- c(sum(wcts!=0)/length(old.vocab), sum(wcts!=0)/new_vocab_size)
  #print some lovely messages
  if(verbose) {
    if(length(new$docs.removed)) {
      cat(sprintf("Removing %i Documents with No Words \n", length(new$docs.removed)))
    }
    
    
    
    cat(sprintf("Your new corpus now has %i documents, %i non-zero terms of %i total terms in the original set. \n%i terms from the new data did not match.\nThis means the new data contained %.1f%% of the old terms\nand the old data contained %.1f%% of the unique terms in the new data. \nYou have retained %i tokens of the %i tokens you started with (%.1f%%).", 
                length(new$documents), #documents
                sum(new$wordcounts>0), #nonzero terms
                length(new$vocab), #total terms
                length(new$words.removed),#unmatched terms
                round(new$prop.overlap[1]*100,3),
                round(new$prop.overlap[2]*100,3),
                sum(new$wordcounts),#total tokens
                total_tokens,#starting tokens
                round(100*sum(new$wordcounts)/total_tokens,3) #percent
    )) 
  }
  new$documents <- lapply(new$documents, function(x) matrix(as.integer(x), nrow=2))
  return(new)
}