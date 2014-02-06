\name{labelTopics}
\alias{labelTopics}
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

"frequency" is the words within each topic with the highest
probability (inferred directly from topic-word distribution parameter \eqn{\beta}).  

"frex" is the words that are both frequent and
exclusive, identifying words that distinguish topics.  This is calculated by taking
the harmonic mean of rank by probability within the topic (frequency) and rank by distribution
of topic given word \eqn{p(z|w=v)} (exclusivity).  In estimating exclusivity we use a 
James-Stein type shrinkage estimator of the distribution \eqn{p(z|w=v)}.

"score" and "lift" are measures provided in two other popular text mining packages. For more
information on type "score", see the R package \code{lda}.  For more
information on type "lift", see Taddy, "Multinomial Inverse Regression
for Text Analysis", Journal of the American Statistical Association 108,
2013 and the R package \code{textir}.
}
\value{
Returns a labelTopics object which includes the words for each label type.}
\seealso{
\code{\link{stm}}
}
\examples{
labelTopics(gadarianFit)
}
