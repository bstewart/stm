\name{plotRemoved}
\alias{plotRemoved}

\title{
Plot documents, words and tokens removed at various word thresholds}
\description{
A plot function which shows the results of using different thresholds in \code{prepDocuments} on the size of the corpus.
}
\usage{
plotRemoved(documents, lower.thresh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{documents}{
    The documents to be used for the stm model
  }
  \item{lower.thresh}{
    A vector of integers, each of which will be tested as a lower threshold 
    for the prepDocuments function.
   }
}
\details{For a lower threshold, \code{prepDocuments} will drop words which appear in fewer than that number of documents, and remove documents which contain no more words. This function allows the user to pass a vector of 
lower thresholds and observe how \code{prepDocuments} will handle each threshold. This function produces  three plots, showing the number of words, the number of documents, and the total number of tokens removed as a function of threshold values. A dashed red line is plotted at the total number of documents, words and tokens respectively.
}
\value{
Invisibly returns a list of
\item{lower.thresh}{The sorted threshold values}
\item{ndocs}{The number of documents dropped for each value of the lower threshold}
\item{nwords}{The number of entries of the vocab dropped for each value of the lower threshold.}
\item{ntokens}{The number of tokens dropped for each value of the lower threshold.}
}
\seealso{
\code{\link{prepDocuments}}
}
\examples{
\dontrun{
plotRemoved(poliblog5k.docs, lower.thresh=seq(from = 10, to = 1000, by = 10))
}
}
