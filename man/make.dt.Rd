\name{make.dt}
\alias{make.dt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a \code{data.table} of topic proportions.
}
\description{
Combines the document-topic loadings (theta) with metadata to create a \code{data.table} object for easy querying.
}
\usage{
make.dt(model, meta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The \code{stm} model.
}
  \item{meta}{
Optionally, the metadata object passed to the \code{stm} model.
}
}
\details{
This is a simple utility function that creates a \pkg{data.table} object which you can use to create 
more complicated queries than via \code{\link{findThoughts}}.  Topics are named via the convention 
\code{Topic#}, for example \code{Topic1}, \code{Topic2} etc.  The object also contains \code{docnum}
which gives the index of the document so you can set keys without worrying about the texts getting
disconnected.
}

\seealso{
\code{\link{findThoughts}}
}
\examples{
dt <- make.dt(gadarianFit, meta=gadarian)
#now we can do any query.  For example the 5 least associated documents with Topic 2 in
#the treated group
dt[treatment==0, docnum[order(Topic2, decreasing=FALSE)][1:5]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
