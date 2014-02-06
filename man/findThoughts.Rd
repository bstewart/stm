\name{findThoughts}
\alias{findThoughts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find Thoughts
}
\description{
  Outputs most representative documents for a particular topic. Use this in order to get a better sense of the content of actual documents with a high topical content.
}
\usage{
findThoughts(model, texts=NULL, topics=NULL, n=3)
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
}

\details{
Prints to the screen the top \code{n} documents ranked by the MAP estimate of the topic's
theta value (which captures the modal estimate of the proportion of word tokens assigned to the topic under the model). Returns document indices and top thoughts.
}
\value{
\item{index}{Matrix with one column per requested topic and one row for
  each requested document.  Entries are the document indices within the dataset.}
\item{docs}{List of top thoughts in the order of topics.}
}
\examples{
findThoughts(gadarianFit, texts=gadarian$open.ended.response, topics=c(1,2), n=3)
}
