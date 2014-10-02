\name{plot.searchK}
\alias{plot.searchK}
\title{
Plots diagnostic values resulting from searchK
}
\description{
 Takes the result of searchK and produces a set of plots for evaluating optimal
 topic numbers via visual representation of diagnostic functions. 
}
\usage{
\method{plot}{searchK}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A searchK object, containing the diagnostic information of an stm with a variety
    of topics.}
\item{...}{
additional arguments for S3 compatability.
}
}
\examples{
\dontrun{
K<-c(5,10,15) 
temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
documents <- out$documents
vocab <- out$vocab
meta <- out$meta
set.seed(02138)
K<-c(5,10,15) 
kresult <- searchK(documents, vocab, K, prevalence=~treatment + s(pid_rep), data=meta)

plot.searchK(kresult)
}
 
}
