#' Find Thoughts
#' 
#' Outputs most representative documents for a particular topic. Use this in
#' order to get a better sense of the content of actual documents with a high
#' topical content.
#' 
#' Returns the top \code{n} documents ranked by the MAP estimate of the topic's
#' theta value (which captures the modal estimate of the proportion of word
#' tokens assigned to the topic under the model). Setting the \code{thresh}
#' argument allows the user to specify a minimal value of theta for returned
#' documents. Returns document indices and top thoughts.
#' 
#' Sometimes you may want to find thoughts which have more conditions than simply 
#' a minimum threshold.  For example, you may want to grab all documents which satisfy
#' certain conditions on the metadata or other topics.  You can supply a query in the
#'  style of \pkg{data.table} to the \code{where} argument.  Note that in \code{data.table}
#'   variables are referenced by their names in the \code{data.table} object.  The topics
#'   themselves are labeled \code{Topic1}, \code{Topic2} etc.  If you supply the metadata
#'    to the \code{meta} argument, you can also query based on any available metadata. 
#'     See below for examples.
#'
#' If you want to pass even more complicated queries, you can use the function \code{\link{make.dt}} 
#' to generate a \code{data.table} object where you can write your own queries.
#'
#' The \code{plot.findThoughts} function is a shortcut for the \code{plotQuote}
#' function.
#' 
#' @aliases findThoughts print.findThoughts plot.findThoughts
#' @param model Model object created by \code{stm}.
#' @param texts A character vector where each entry contains the text of a
#' document.  Must be in the same order as the documents object. NOTE: This is not the
#' documents which are passed to \code{stm} and come out of \code{prepDocuments}, 
#' this is the actual text of the document.
#' @param topics The topic number or vector of topic numbers for which you want
#' to find thoughts.  Defaults to all topics.
#' @param n The number of desired documents to be displayed per topic.
#' @param thresh Sets a minimum threshold for the estimated topic proportion
#' for displayed documents.  It defaults to imposing no restrictions.
#' @param where An expression in the form of a \code{data.table} query. This is passed to the \code{i} argument in data.table and a custom query is passed to \code{j}.  This cannot be used with \code{thresh}.  See below for more details.
#' @param meta The meta data object to be used with \code{where}.
#' @return A \code{findThoughts} object:
#' \item{index}{List with one entry per
#' topic.  Each entry is a vector of document indices.} 
#' \item{docs}{List with
#' one entry per topic.  Each entry is a character vector of the corresponding
#' texts.}
#' @seealso \code{\link{plotQuote}}
#' @examples
#' findThoughts(gadarianFit, texts=gadarian$open.ended.response, topics=c(1,2), n=3)
#' 
#' #We can plot findThoughts objects using plot() or plotQuote
#' thought <- findThoughts(gadarianFit, texts=gadarian$open.ended.response, topics=1, n=3)
#' 
#' #plotQuote takes a set of sentences
#' plotQuote(thought$docs[[1]])
#' 
#' #we can use the generic plot as a shorthand which will make one plot per topic
#' plot(thought)
#' 
#' #we can select a subset of examples as well using either approach
#' plot(thought,2:3)
#' plotQuote(thought$docs[[1]][2:3])
#' 
#' 
#' #gather thoughts for only treated documents
#' thought <- findThoughts(gadarianFit, texts=gadarian$open.ended.response, topics=c(1,2), n=3, 
#'                        where = treatment==1, meta=gadarian)
#' plot(thought)
#' #you can also query in terms of other topics
#' thought <- findThoughts(gadarianFit, texts=gadarian$open.ended.response, topics=c(1), n=3, 
#'                         where = treatment==1 & Topic2>.2, meta=gadarian)
#' plot(thought)         
#' #these queries can be really complex if you like
#' thought <- findThoughts(gadarianFit, texts=gadarian$open.ended.response, topics=c(1), n=3, 
#'                        where = (treatment==1 | pid_rep > .5) & Topic3>.2, meta=gadarian)
#' plot(thought)         
#' @export
#' @import data.table
findThoughts <- function(model, texts=NULL, topics=NULL, n=3, thresh=NULL,
                         where=NULL, meta=NULL) {
#Grab up to n texts which are above the threshold
  theta <- model$theta
  if(is.null(topics)) topics <- 1:ncol(theta)
  if(!is.null(texts) && length(texts)!=nrow(theta)) stop("Number of provided texts and number of documents modeled do not match")
  if(!is.null(texts) && class(texts)=="list" && class(texts[[1]])=="matrix") stop("It looks like you are trying to pass the numeric documents object. \n The texts argument wants a vector of characters that contain the actual text of the documents.")
  #since we are allowing n to be infinity have to recode if its too high.
  if(n > nrow(theta)) n <- nrow(theta)
  if(n < 1) stop("Must request at least one returned document.")
  out <- list()
  where <- substitute(where)

    #First check if they are trying to use thresh and where together
  if(!is.null(where) & !is.null(thresh)) warning("Threshold value is ignored when where argument is non null \n Include threshold explicitly in where statement.")
  #okay now if thresh isn't used fill it in as 0
  if(is.null(thresh)) thresh <- 0.0
  
  #If they are using "where" we need the data.table infrastructure
  if(!is.null(where)) {
    colnames(theta) <- sprintf("Topic%i", 1:ncol(theta))
    if(is.null(meta)) {
      dt <- data.table::data.table(docnum=1:nrow(theta),theta)
    } else {
      dt <- data.table::data.table(docnum=1:nrow(theta),theta, meta)
    }
    for(i in 1:length(topics)) {
      what <- parse(text=sprintf("docnum[order(Topic%i, decreasing=TRUE)][1:%i]", topics[i], n))
      index <- dt[eval(where), eval(what)]
      out$index[[i]] <- index
      #grab the associated texts
      if(!is.null(texts)) out$docs[[i]] <- texts[index]
    }
  } else {
    #Here we are in the simpler non-"where" branch of the tree
    for(i in 1:length(topics)) {
      k <- topics[i]
      #grab the values and the rank
      index <- order(theta[,k], decreasing=TRUE)[1:n]
      val <- sort(theta[,k], decreasing=TRUE)[1:n]
      #subset to those values which meet the threshold
      index <- index[which(val>=thresh)]
      out$index[[i]] <- index
      #grab the associated texts
      if(!is.null(texts)) out$docs[[i]] <- texts[index]
    }
  }
  names(out$index) <- paste("Topic", topics)
  if(!is.null(texts)) names(out$docs) <- paste("Topic", topics)
  
  class(out) <- "findThoughts"
  
  if(is.null(texts)) return(out)
  
  out$docs <- lapply(out$docs, unlist)
  if(is.factor(texts)){
    warning("texts are of type 'factor.'  Converting to character vectors.  Use 'as.character' to avoid this warning in the future.")
    out$docs <- lapply(out$docs, as.character)  
  }
  
  return(out)
}

#' @method print findThoughts
#' @export
print.findThoughts <- function(x,...) {
  toprint <- vector(length=length(x$docs))
  for(i in 1:length(x$docs)) {
    docs <- paste(x$docs[[i]], collapse="\n \t")
    toprint[i] <- sprintf("\n %s: \n \t %s", names(x[[1]])[i], docs)
  }
  cat(toprint)
}

#' @method plot findThoughts
#' @export
plot.findThoughts <- function(x, sentences=NULL, ...) {
  for(i in 1:length(x$docs)) {
    if(is.null(sentences)) sentences <- 1:length(x$docs[[i]])
    plotQuote(x$docs[[i]][sentences],...)
  }
}

#' Make a \code{data.table} of topic proportions.
#' 
#' Combines the document-topic loadings (theta) with metadata to create a \code{data.table} object for easy querying.
#' 
#' This is a simple utility function that creates a \pkg{data.table} object which you can use to create 
#' more complicated queries than via \code{\link{findThoughts}}.  Topics are named via the convention 
#' \code{Topic#}, for example \code{Topic1}, \code{Topic2} etc.  The object also contains \code{docnum}
#' which gives the index of the document so you can set keys without worrying about the texts getting
#' disconnected.
#' 
#' We expect that for the vast majority of users the functionality in \code{\link{findThoughts}} will be
#' sufficient.
#' 
#' @param model The \code{stm} model.
#' @param meta Optionally, the metadata object passed to the \code{stm} model.
#' @seealso \code{\link{findThoughts}}
#' @examples
#' dt <- make.dt(gadarianFit, meta=gadarian)
#' #now we can do any query.  For example the 5 least associated documents with Topic 2 in
#' #the treated group
#' dt[treatment==0, docnum[order(Topic2, decreasing=FALSE)][1:5]]
#' 
#' @export
make.dt <- function(model, meta=NULL) {
  theta <- model$theta
  colnames(theta) <- sprintf("Topic%i", 1:ncol(theta))
  if(is.null(meta)) {
    return(data.table::data.table(docnum=1:nrow(theta),theta))
  } else {
    return(data.table::data.table(docnum=1:nrow(theta),theta, meta))
  }
}


