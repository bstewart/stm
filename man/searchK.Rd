\name{searchK}
\alias{searchK}
\title{
Computes diagnostic values for models with different values of K (number of topics).}

\description{
 With user-specified initialization, this function  runs selectModel for different
 user-specified topic numbers and computes diagnostic properties for the resturned model. These include exclusivity,
 semantic coherence, heldout likelihood, bound, lbound, and residual. 
}
\usage{
 searchK(documents, vocab, K, init.type = "Spectral", 
                    N=floor(.1*length(documents)), proportion=.5, heldout.seed=NULL,
                    M=10,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{documents}{
    The documents to be used for the stm model
  }
  \item{vocab}{
    The vocabulary to be used for the  stmmodel
    }
   \item{K}{
    A vector of different topic numbers
}
\item{init.type}{
    The method of initialization.  Must be either Latent
    Dirichlet Allocation (LDA), Dirichlet Multinomial Regression
    Topic Model (DMR), a random initialization, `spectral', or a previous STM
    object. If an initialization other than spectral is used, then the selection procedure is manyTopics is used. 
}
\item{N}{
  Number of docs to be partially held out
}
\item{proportion}{
  Proportion of docs to be held out.
}   

\item{heldout.seed}{
  If desired, a seed to use when holding out documents for later heldout likelihood
  computation
} 

\item{M}{
  M value for exclusivity computation
}

\item{...}{
    Other diagnostics parameters.
  }


}

\value{
\item{exclus}{Exclusivity of each model.}
\item{semcoh}{Semantic coherence of each model.}
\item{heldout}{Heldout likelihood for each model.}
\item{residual}{Residual for each model.}
\item{bound}{Bound for each model.}
\item{lbound}{lbound for each model.}
\item{em.its}{Total number of EM iterations used in fiting the model.}
}

\examples{


\dontrun{

K<-c(5,10,15) 
temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
documents <- out$documents
vocab <- out$vocab
meta <- out$meta
set.seed(02138)
K<-c(5,10,15) 
kresult <- searchK(documents, vocab, K, prevalence=~treatment + s(pid_rep), data=meta)


}
 
}
