\name{plotRemoved}
\alias{plotRemoved}

\title{
Produces a plot of number of words, documents, and tokens removed by prepDocuments for different threshholds.}
\description{
For a lower threshold, prepDocuments will drop words which appear in fewer than that number of documents,
and remove documents which contain no more words. countRemoved allows the user to pass a vector of 
lower thresholds and observe how prepDocuments will handle each threshold. This function produces 
three plots, showing the number of words, the number of documents, and the total number of tokens
removed as a function of threshold values. 
}
\usage{
plotRemoved(documents, vocab, meta, lower.thresh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{documents}{
    The documents to be used for the stm model
  }
  \item{vocab}{
    The vocabulary to be used for the stm model
    }
   \item{meta}{
    The metadata to be used for the stm model
   }
   \item{lower.thresh}{
    A vector of integers, each of which will be tested as a lower threshold 
    for the prepDocuments function.
   }
}

\details{
This function produces a plot, and does not return values.
}

\examples{
\dontrun{
lower.thresh<-seq(from = 10, to = 1000, by = 10)
plotRemoved(docs, vocab, meta, lower.thresh)
}
}
