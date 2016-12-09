#rewritten from scratch on 11/7/2014
#new version is dramatically faster but still not as memory efficient
#as it could be.  I'm restoring and recalculating the thresholds at every
#level but theoretically this could all be constructed as a cumulative
#measure.  Still it would take an enormous document set for this to matter.
# -BMS
plotRemoved<- function(documents, lower.thresh){
  #sort so we know it is in order
  lower.thresh <- sort(lower.thresh)
  
  ###
  #Create some useful representations
  ###
  #standard triplet form
	triplet <- doc.to.ijv(documents)
  #calculate the number of docs for each vocab item
  wordcounts <- tabulate(triplet$j)
  #calculate the number of tokens for each vocab item
  tokencount <- tabulate(rep(triplet$j, times=triplet$v))              
  
  ###
  # Calculate quantities of interest
  ###
  #which words will drop?
  drop <- sapply(lower.thresh, function(x) which(wordcounts<=x), simplify=FALSE)
  
  #number of words dropped is just the length
  nwords <- unlist(lapply(drop, length))
  #for tokens we sum over the token counts for the dropped words
  ntokens <- unlist(lapply(drop, function(x) sum(tokencount[x])))
  
  #for documents its a bit more nuanced...
  
  #calculate the number of docs in which the word with the highest count
  #appears.  this tells us what number we would have to drop to lose the doc.
  docthresh <- unlist(lapply(documents, function(x) max(wordcounts[x[1,]])))
  ndocs <- sapply(lower.thresh, function(x) sum(docthresh<=x), simplify=TRUE)
  
	# Composite Plot
  oldpar <- par(no.readonly=TRUE)
	par(mfrow = c(1,3), oma = c(2,2,2,2))
	plot(lower.thresh, ndocs, type = "n", xlab = "", 
       ylab = "Number of Documents Removed", main = "Documents Removed by Threshold")
	lines(lower.thresh, ndocs, lty=1, col=1) 
	abline(a = length(documents), lty=2, b=0, col="red")

	plot(lower.thresh, nwords, type = "n", 
       xlab = "Threshold (Minimum No. Documents Appearing)", 
       ylab = "Number of Words Removed", main = "Words Removed by Threshold")
	lines(lower.thresh, nwords, lty=1, col=1) 
	abline(a = length(tokencount), lty=2, b=0, col="red")
  
  plot(lower.thresh, ntokens, type = "n", xlab= "", 
       ylab = "Number of Tokens Removed", main = "Tokens Removed by Threshold")
  lines(lower.thresh, ntokens, lty=1, col=1)
  abline(a = sum(tokencount), lty=2, b=0, col="red")
  par(oldpar)
  return(invisible(list(lower.thresh=lower.thresh,
                        ndocs=ndocs,
                        nwords=nwords,
                        ntokens=ntokens)))
}