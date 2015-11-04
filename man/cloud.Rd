\name{cloud}
\alias{cloud}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a wordcloud
}
\description{
Use the \pkg{wordcloud} package to plot a wordcloud for a particular topic
}
\usage{
cloud(stmobj, topic=NULL, type=c("model", "documents"), 
      documents, thresh=.9, max.words=100, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stmobj}{
The STM model object to be used in making the word cloud.
}
\item{topic}{
NULL to plot the marginal distribution of words in the corpus, or a single integer indicating the topic number.
}
\item{type}{
Specifies how the wordcloud is constructed.  The type \code{"model"} which is used by default is based on the probability of the word given the topic.  The type \code{"documents"} plots words within documents that have a topic proportion of higher than \code{thresh}.  This requires that the \code{documents} argument also be specified.   
}
\item{documents}{
The documents object of the same kind as passed to \code{\link{stm}}.  This is only necessary if \code{type="documents"}.
}
\item{thresh}{
The threshold for including a document in the \code{type="documents"} setting.
}
\item{max.words}{
The maximum number of words to be plotted.
}
\item{...}{
Additional parameters passed to \code{wordcloud}.
}
}
\details{
Uses the \pkg{wordcloud} package to make a word cloud of a particular topic.  The option \code{"model"} uses the topic-word model parameters.  Thus it shows words weighted by their probability conditional that the word comes from a particular topic.  With content covariates it averages over the values for all levels of the content covariate weighted by the empirical frequency in the dataset.  The option \code{"documents"} plots the words which appear in documents that have a topic proportion higher than \code{thresh}.  Thus \code{"model"} gives a pure model based interpetation of the topic while \code{"documents"} gives a picture of all the words in documents which are highly associated with the topic.
}
\references{
Ian Fellows (2014). wordcloud: Word Clouds. R package version 2.5.
  \url{http://CRAN.R-project.org/package=wordcloud}}
\seealso{
\code{\link{plot.STM}}
}
\examples{
cloud(gadarianFit, 1)
}
