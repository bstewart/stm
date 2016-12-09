\name{gadarian}
\alias{gadarian}
\alias{gadarianFit}
\docType{data}
\title{
Gadarian and Albertson data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This data set contains variables from Gadarian and Albertson (2014).
The experiment had those in the treatment condition write about what made them anxious about immigration. The control condition just had subjects write about immigration. 
}
\usage{data(gadarian)}
\format{
  A data frame with 351 observations on the following 3 variables.
  \describe{
    \item{\code{MetaID}}{a numeric vector containing identification numbers; not used for analysis}
    \item{\code{treatment}}{a numeric vector indicating treatment condition}
    \item{\code{pid_rep}}{a numeric vector of party identification}
    \item{\code{open.ended.response}}{a character vector of the subject's open ended response}
  }
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Gadarian, Shana Kushner, and Bethany Albertson. "Anxiety, immigration, and the search for information." Political Psychology 35.2 (2014): 133-164.

Roberts, Margaret E., Brandon M. Stewart, Dustin Tingley, Christopher Lucas, Jetson Leder-Luis, Shana Kushner Gadarian, Bethany Albertson, and David G. Rand.  "Structural Topic Models for Open-Ended Survey Responses." American Journal of Political Science 58, no 4 (2014): 1064-1082.

}
\examples{
head(gadarian)
#Process the data for analysis.
temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
meta<-temp$meta
vocab<-temp$vocab
docs<-temp$documents
out <- prepDocuments(docs, vocab, meta)
docs<-out$documents
vocab<-out$vocab
meta <-out$meta
}
\keyword{datasets}
