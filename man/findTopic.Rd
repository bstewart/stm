\name{findTopic}
\alias{findTopic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find topics that contain user specified words.
}
\description{
Find topics that contain user specified words.
}
\usage{
findTopic(x, list, n=20, type=c("prob", "frex", "lift","score"), verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The STM model object to be searched. May also be the output from sageLabels.
}
\item{list}{
Character vector containing words to be searched.
}
\item{n}{
number of words to consider
}
\item{type}{
Type of words to be searched.
}
\item{verbose}{
a logical indicating whether details should be printed to the screen.
}
}
\examples{
\dontrun{
lab <- sageLabels(gadarianFit, 5)
findTopic(lab, c("poor", "immigr", "peopl"))
findTopic(gadarianFit, c("poor", "immigr", "peopl"))
}
}
