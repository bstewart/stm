\name{gadarian}
\alias{gadarian}
\alias{gadarianFit}
\docType{data}
\title{
Gadarian and Albertson data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This data set contains variables from Gadarian and Albertson, forthcoming, "Anxiety, Immigration, and the Search for Information", Political Psychology.
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
Gadarian and Albertson, forthcoming, "Anxiety, Immigration, and the Search for Information", Political Psychology
http://onlinelibrary.wiley.com/doi/10.1111/pops.12034/abstract

Roberts, M., Stewart, B., Tingley, D., Lucas, C., Leder-Luis, J., Gadarian, S., Albertson, B., Albertson, B. and Rand, D. (Forthcoming). "Structural topic models for open ended survey responses." American Journal of Political Science
http://scholar.harvard.edu/files/dtingley/files/topicmodelsopenendedexperiments.pdf
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
