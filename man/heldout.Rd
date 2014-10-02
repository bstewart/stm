\name{make.heldout}
\alias{make.heldout}
\alias{eval.heldout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heldout Likelihood by Document Completion}
\description{
Tools for making and evaluating heldout datasets.
}
\usage{
make.heldout(documents, vocab, 
             N=floor(.1*length(documents)), 
             proportion=.5, seed=NULL)
eval.heldout(model, missing)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{documents}{the documents to be modeled.
}
  \item{vocab}{
  the vocabulary item
  }
  \item{N}{
  number of docs to be partially held out
  }
  \item{proportion}{
  proportion of docs to be held out.
  }
  \item{seed}{
  the seed, set for replicability
  }
  \item{model}{
  an stm model
  }
  \item{missing}{
  a missing object created by make.heldout
  }
}
\details{
These functions are used to create and evaluate heldout likelihood using the document completion method.  The basic idea is to hold out some fraction of the words in a set of documents, train the model and use the document-level latent variables to evaluate the probability of the heldout portion. See the example for the basic workflow.  
}
\examples{
\dontrun{
prep <- prepDocuments(poliblog5k.docs, poliblog5k.voc, 
                      poliblog5k.meta,subsample=500,
                      lower.thresh=20,upper.thresh=200)
heldout <- make.heldout(prep$documents, prep$vocab)
documents <- heldout$documents
vocab <- heldout$vocab
meta <- out$meta

stm1<- stm(documents, vocab, 5, 
           prevalence =~ rating+ s(day), 
           init.type="Random",
           data=meta, max.em.its=5)
eval.heldout(stm1, heldout$missing)
}
}
