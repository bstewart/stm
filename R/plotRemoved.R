plotRemoved<- function(documents, vocab, meta, lower.thresh){
	# Initialize vector of removed words to be lenght of lower.thresh
	wordsremoved <- vector("numeric", length(lower.thresh))
	docsremoved <- vector("numeric", length(lower.thresh))
  tokensremoved <- vector("numeric", length(lower.thresh))

	# Loop over lower.thresh values
	for(i in 1:length(lower.thresh)){
		lt <- lower.thresh[i] # Grab lower.thresh value

		#Run prep-docs with that lower threshold
		temp<- prepDocuments(documents, vocab, meta, lower.thresh = lt, verbose = FALSE)

		#Store number of words removed
		docsremoved[i] <- length(temp$docs.removed)
		wordsremoved[i] <- length(temp$words.removed)
    tokensremoved[i] <- temp$tokens.removed
	}

	# Snip when reaches max
	wordsstop <- match(length(vocab), wordsremoved)
	docsstop <- match(length(documents), docsremoved)
  tokensstop <- match(sum(temp$wordcounts), tokensremoved)

	if(!is.na(wordsstop)){
		wordsremoved <- wordsremoved[1:wordsstop]
		wordsthresh <- lower.thresh[1:wordsstop]
	}else{
		wordsthresh <- lower.thresh
	}
	if(!is.na(docsstop)){
		docsremoved <- docsremoved[1:docsstop]
		docsthresh <- lower.thresh[1:docsstop]
	}else{
		docsthresh <- lower.thresh
	}
  if(!is.na(tokensstop)){
    tokensremoved <- tokensremoved[1:tokensstop]
    tokensthresh <- lower.thresh[1:tokensstop]
  }else{
    tokensthresh <- lower.thresh
  }

	# Composite Plot
	par(mfrow = c(1,3), oma = c(2,2,2,2))
	plot(docsthresh, docsremoved, type = "n", xlab = "", 
       ylab = "Number of Documents Removed", main = "Docs Removed by Threshold")
	lines(docsthresh, docsremoved, lty=1, col=1) 
	abline(a = length(documents), lty=2, b=0, col="red")

	plot(wordsthresh, wordsremoved, type = "n", xlab = "Threshold (Minimum No. Documents Appearing)", 
       ylab = "Number of Words Removed", main = "Words Removed by Threshold")
	lines(wordsthresh, wordsremoved, lty=1, col=1) 
	abline(a = length(vocab), lty=2, b=0, col="red")
  
  plot(tokensthresh, tokensremoved, type = "n", xlab= "", 
       ylab = "Number of Tokens Removed", main = "Tokens Removed by Threshold")
  lines(tokensthresh, tokensremoved, lty=1, col=1)
  abline(a = sum(temp$wordcounts), lty=2, b=0, col="red")
}