\name{termitewriter}
\alias{termitewriter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Termite writer
}
\description{
Writes output necessary to use the termite d3 visualization.}
\usage{
termitewriter(app.path, model, meta, DocID=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{app.path}{
    Path to the output folder, which termite can then use to load the documents
    }
  \item{model}{
    STM model output object.
    }
   \item{meta}{
     Metadata with covariates the user wants to vizualize.
   }
   \item{DocID}{
     If the user wants to associate each document with a particular
     docID, they can input that here.  By default, the program will
     create a sequential document ID.
     }
}
\examples{
#termitewriter("~/Out", mod.out, data) 
}

