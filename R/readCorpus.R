

###
# Corpus Readers
###

readCorpus <- function(corpus, type=c("dtm", "ldac", "slam", "Matrix","txtorgvocab")) {
  type <- match.arg(type)
  switch(type,
         dtm = read.dtm(corpus),
         ldac = read.ldac(corpus),
         slam = read.slam(corpus),
         txtorgvocab = read.txtorg.vocab(corpus),
         Matrix = read.slam(as.simple_triplet_matrix(corpus)))
}

read.ldac <- function(filename) {
  #Read the .ldac format
  # Based on Jonathan Chang's  reader with addition of zero correction.
  d <- scan(filename, what = "character", sep = "\n")
  d <- chartr(":", " ", d)
  d <- strsplit(d, " ", fixed = TRUE)
  d <- lapply(d, function(x) matrix(as.integer(x[-1]), nrow = 2))
  mapply(function(x) rbind(x[1,]+1, x[2,]), d) #zero correction
}

read.dtm <- function(dtm) {
  #test for and adjust for mispecification
  if("simple_triplet_matrix" %in% class(dtm)) {
    warning("Please use the slam option.  dtm is for dense matrices.")
    read.slam(dtm)
  }
  #convert a standard document-term matrix to list format.
  dtm.mat <- as.matrix(dtm)
  vocab <- colnames(dtm)
  documents <- apply(dtm.mat, 1, function(y) {
    rbind(which(y > 0), as.integer(y[y > 0])) })
  return(list(documents=documents, vocab=vocab))
}

read.slam <- function(corpus) {
  #convert a simple triplet matrix to list format.
  if(!inherits(corpus, "simple_triplet_matrix")) stop("corpus is not a simple triplet matrix")
  documents <- ijv.to.doc(corpus$i, corpus$j, corpus$v)
  vocab <- corpus$dimnames[[2]]
  return(list(documents=documents,vocab=vocab))
}

read.txtorg.vocab <- function(filename) {
  return(readLines(filename, encoding="UTF-8",warn=FALSE))
}