convertCorpus <- function(documents, vocab, type=c("slam", "lda", "Matrix")) {
  #a function to output our format to one of the existing frameworks
  type <- match.arg(type)
  switch(type,
         slam = convert.slam(documents,vocab),
         lda = convert.lda(documents,vocab),
         Matrix = convert.Matrix(documents, vocab))
}

convert.slam <- function(documents, vocab) {
  slam <- doc.to.ijv(documents)
  slam <- simple_triplet_matrix(slam$i, slam$j, slam$v)
  colnames(slam) <- vocab
  return(slam)
}

convert.lda <- function(documents,vocab) {
  #just subtract one off from the indices
  documents<- lapply(documents, function(x) {
    x[1, ] <- as.integer(x[1, ] - 1)
    x})
  return(list(documents=documents, vocab=vocab))
}

convert.Matrix <- function(documents,vocab) {
  docs <- doc.to.ijv(documents)
  mat <- sparseMatrix(docs$i,docs$j, x=docs$v)
  colnames(mat) <- vocab
  return(mat)
}

