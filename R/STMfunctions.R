####
# Non-Exported Utility Functions
####

##
# Random Utilities

#function for collapsing a character vector to a comma separated list
#note that we eliminate things with zero length so that we can return
#fewer than n words and still have the lists look nice
commas <- function(text){  
  paste(text[nchar(text)>0], collapse=", ")
}

#from the R documentation for is.integer
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
} 

posint <- function(x) {
  all(is.wholenumber(x)) & all(x>0)
}

#From Peter Hoff
rmvnorm<-function(n,mu,Sigma,chol.Sigma=chol(Sigma)) {
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol.Sigma) +c(mu))
}

##
# Label Functions
##
calcfrex <- function(logbeta, w=.5, wordcounts=NULL) {
  excl <- t(t(logbeta) - col.lse(logbeta))
  if(!is.null(wordcounts)) {
    #if word counts provided calculate the shrinkage estimator
    excl <- safelog(sapply(1:ncol(excl), function(x) js.estimate(exp(excl[,x]), wordcounts[x])))
  } 
  freqscore <- apply(logbeta,1,rank)/ncol(logbeta)
  exclscore <- apply(excl,1,rank)/ncol(logbeta)
  frex <- 1/(w/freqscore + (1-w)/exclscore)
  apply(frex,2,order,decreasing=TRUE)
}

#A James-Stein Estimator Shrinking to a Uniform Distribution
#This draws from the Hausser and Strimmer (2009) JMLR piece.
js.estimate <- function(prob, ct) {
  # MLE of prob estimate
  mlvar <- prob*(1-prob)/(ct-1)
  unif <- 1/length(prob) 
  
  # Deviation from uniform
  deviation <- sum((prob-unif)^2)
  
  #take care of special case,if no difference it doesn't matter
  if(deviation==0) return(prob)
  
  lambda <- sum(mlvar)/deviation
  #truncate
  if(lambda>1) lambda <- 1
  if(lambda<0) lambda <- 0
  
  #Construct shrinkage estimator as convex combination of the two
  lambda*unif + (1 - lambda)*prob
}

#Calculate Lift (Taddys thing this is beta_k,v divided by the empirical term probability)
calclift <- function(logbeta, wordcounts) {
  emp.prob <- log(wordcounts) - log(sum(wordcounts))
  lift <- logbeta - rep(emp.prob, each=nrow(logbeta)) 
  apply(lift, 1, order, decreasing=TRUE)
}

#Calculate score (Chang in LDA package etc.)
#beta times (logbeta - mean_k(logbeta_v))
calcscore <- function(logbeta) { 
  ldascore <- exp(logbeta)*(logbeta - rep(colMeans(logbeta), each=nrow(logbeta)))
  apply(ldascore, 1, order, decreasing=TRUE)
} 
##
#Document convertors
##

#Our Format to Triplet Format
doc.to.ijv <- function(documents, fixzeroindex=TRUE) {
  #Turns our format into triplet format (can be zero indexed)
  indices <- unlist(lapply(documents, '[',1,)) #grab the first row
  if((0 %in% indices) & fixzeroindex) indices <- indices + 1 #if zero-indexed, fix it.
  counts <- lapply(documents, '[',2,)  #grab the second row but preserve the list structure for a moment
  VsubD <- unlist(lapply(counts,length)) #grab the number of unique words per document
  rowsums <- unlist(lapply(counts,sum)) #grab the number of tokens per documents
  docids <- rep(1:length(documents), times=VsubD) #add row numbers
  counts <- unlist(counts) #unlist the count structure
  #now we return using the convention for the simple triplet matrix,
  #plus the row sums which we use in DMR.
  return(list(i=as.integer(docids), j=as.integer(indices), v=as.integer(counts), rowsums=as.integer(rowsums)))
}

#Triplet Format to our Document Format
ijv.to.doc <- function(i,j,v) {
  index <- split(j,i)
  index <- lapply(index,as.integer)
  count <- split(v,i)
  count <- lapply(count,as.integer)
  mapply(rbind,index,count)
}

##
#A series of fast softmax functions mostly wrappers around matrixStats package functions
##
logsoftmax <- function(x) {
  x - lse(x)
}

lse <- function(x) {
 logSumExp(x)
}

#Note here we have simply cut the NA checks, the names element etc.
row.lse <- function(mat) {
  .Call("rowLogSumExps", mat, FALSE, FALSE, TRUE, PACKAGE = "matrixStats")
  #rowLogSumExps(mat)
}
col.lse <- function(mat) {
  #colLogSumExps(mat)
  .Call("rowLogSumExps", mat, FALSE, FALSE, FALSE, PACKAGE = "matrixStats")
}

softmax <- function(x) {
  exp(x - lse(x))
}

safelog <- function(x) {
  out <- log(x)
  out[which(out< -1000)] <- -1000
  out
}

