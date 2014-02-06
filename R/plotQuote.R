plotQuote <- function(sentences, width=30) {
  xaxt <- yaxt <-"n"
  xlab <- ylab <- ""
  xlim <- c(0,5)
  ylim <- NULL
  numlines=NULL
  out <- list()
  for(j in 1:length(sentences)){
    sentence <- sentences[j]
    out[[j]] <- str_wrap(sentence, width)
    numlines[j] <- length(strsplit(out[[j]], "\n")[[1]])
  }
  if(is.null(ylim)) ylim=c(0,.5*sum(numlines))
  plot(c(0,0), col="white", xlim=xlim, ylim=ylim, xaxt=xaxt, yaxt=yaxt, xlab=xlab, ylab=ylab)
  numlines <- c(0, rev(numlines))
  out <- rev(out)
  start <- 0
  for(j in 2:(length(out)+1)){
    text(2.5, start + numlines[j-1]/4 + numlines[j]/4, out[[j-1]])
    start <- start + numlines[j-1]/4 + numlines[j]/4 
    if(j!=1 & j!=(length(out)+1)) lines(c(0,5), c(sum(numlines[1:j-1])/2 + numlines[j]/2, sum(numlines[1:j-1])/2 + numlines[j]/2), lty=2)
  }
}
