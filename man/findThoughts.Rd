\name{findThoughts}
\alias{findThoughts}
\alias{print.findThoughts}
\alias{plot.findThoughts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find Thoughts
}
\description{
  Outputs most representative documents for a particular topic. Use this in order to get a better sense of the content of actual documents with a high topical content.
}
\usage{
findThoughts(model, texts=NULL, topics=NULL, n=3, thresh=0.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    Model object created by \code{stm}.
    }
    \item{texts}{
    A character vector where each entry contains the text of a document.  Must be in the same order as the documents object.
  }
  \item{topics}{
     The topic number or vector of topic numbers for which you want to find thoughts.  Defaults to all topics.
}
  \item{n}{
    The number of desired documents to be displayed per topic. 
  }
  \item{thresh}{
   Sets a minimum threshold for the estimated topic proportion for displayed documents.  It defaults to
   imposing no restrictions.
  }
}

\details{
Returns the top \code{n} documents ranked by the MAP estimate of the topic's
theta value (which captures the modal estimate of the proportion of word tokens assigned to the topic under the model). Setting the \code{thresh} argument allows the user to specify a minimal value of theta for returned documents. Returns document indices and top thoughts.

The \code{plot.findThoughts} function is a shortcut for the \code{plotQuote} function.
}
\value{
A \code{findThoughts} object
\item{index}{List with one entry per topic.  Each entry is a vector of document indices.}
\item{docs}{List with one entry per topic.  Each entry is a character vector of the corresponding texts.}
}
\examples{
findThoughts(gadarianFit, texts=gadarian$open.ended.response, topics=c(1,2), n=3)

#We can plot findThoughts objects using plot() or plotQuote
thought <- findThoughts(gadarianFit, texts=gadarian$open.ended.response, topics=1, n=3)

#plotQuote takes a set of sentences
plotQuote(thought$docs[[1]])

#we can use the generic plot as a shorthand which will make one plot per topic
plot(thought)

#we can select a subset of examples as well using either approach
plot(thought,2:3)
plotQuote(thought$docs[[1]][2:3])

}
\seealso{
\code{\link{plotQuote}}
}
