\name{labelTopics}
\alias{labelTopics}
\alias{print.labelTopics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Label topics
}
\description{
Generate a set of words describing each topic from a fitted STM object.  Uses a variety of labeling algorithms (see details).
}
\usage{
labelTopics(model, topics=NULL, n = 7, frexweight = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
 An \code{STM} model object.
}
\item{topics}{
  A vector of numbers indicating the topics to include.  Default is all topics.
  }
  \item{n}{
The desired number of words (per type) used to label each topic.
}
  \item{frexweight}{
A weight used in our approximate FREX scoring algorithm (see details).}
}
\details{
Four different types of word weightings are printed with label topics.

Highest Prob: are the words within each topic with the highest
probability (inferred directly from topic-word distribution parameter \eqn{\beta}).  

FREX: are the words that are both frequent and
exclusive, identifying words that distinguish topics.  This is calculated by taking
the harmonic mean of rank by probability within the topic (frequency) and rank by distribution
of topic given word \eqn{p(z|w=v)} (exclusivity).  In estimating exclusivity we use a 
James-Stein type shrinkage estimator of the distribution \eqn{p(z|w=v)}.

Score and Lift are measures provided in two other popular text mining packages. For more
information on type Score, see the R package \code{\link{lda}}.  For more
information on type Lift, see Taddy, "Multinomial Inverse Regression
for Text Analysis", Journal of the American Statistical Association 108,
2013 and the R package \code{textir}.
}
\value{
A labelTopics object (list)
\item{prob }{matrix of highest probability words}
\item{frex }{matrix of highest ranking frex words}
\item{lift }{matrix of highest scoring words by lift}
\item{score }{matrix of best words by score}
\item{topicnums }{a vector of topic numbers which correspond to the rows}
}
\references{
Taddy, Matt. "Multinomial inverse regression for text analysis." Journal of the American Statistical Association 108.503 (2013): 755-770.
}
\seealso{
\code{\link{stm}}
\code{\link{plot.STM}}
}
\examples{
labelTopics(gadarianFit)
}
