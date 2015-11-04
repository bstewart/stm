\name{convertCorpus}
\alias{convertCorpus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert \pkg{stm} formatted documents to another format
}
\description{
Takes an \pkg{stm} formatted documents and vocab object and returns formats useable in other packages.
}
\usage{
convertCorpus(documents, vocab, type=c("slam", "lda", "Matrix"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{documents}{
the documents object in \pkg{stm} format
}
\item{vocab}{
the vocab object in \pkg{stm} format
}
\item{type}{
the output type desired.  See Details.
}
}
\details{
The various type conversions are described below:
  \describe{
    \item{\code{type="slam"}}{Converts to the simple triplet matrix representation used by the \pkg{slam} package.  This is the format used internally by \pkg{tm}.}
    \item{\code{type="lda"}}{Converts to the format used by the \pkg{lda} package.  This is a very minor change as the format in \pkg{stm} is based on \pkg{lda}'s data representation.  The difference as noted in \code{\link{stm}} involves how the numbers are indexed.  Accordingly this type returns a list containing the new documents object and the unchanged vocab object.}
    \item{\code{type="Matrix"}}{Converts to the sparse matrix representation used by \pkg{Matrix}.  This is the format used internally by numerous other text analysis packages.}
}
If you want to write out a file containing the sparse matrix representation popularized by David Blei's \code{C} code \code{ldac} see the function \code{\link{writeLdac}}.
}

\seealso{
\code{\link{writeLdac}}
\code{\link{readCorpus}}
\code{\link{poliblog5k}}
}
\examples{
#convert the poliblog5k data to slam package format
poliSlam <- convertCorpus(poliblog5k.docs, poliblog5k.voc, type="slam")
class(poliSlam)
}
