\name{sageLabels}
\alias{sageLabels}
\alias{print.sageLabels}
\title{
Displays verbose labels that describe topics and topic-covariate groups in depth.}
\description{
For each topic or, when there is a covariate at the bottom of the model, for each topic-covariate
group, sageLabels provides a list of the highest marginal probability words, the highest marginal
FREX words, the highest marginal lift words, and the highest marginal
score words, where marginal means it is summing over all potential covariates.  It also 
provides each topic's Kappa (words associated with each topic) and
baselined Kappa (baseline word distribution).
}
\usage{
sageLabels(model, n=7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
       \item{model}{
      A fitted STM model object.
      }

      \item{n}{
	     The number of words to print per topic/topic-covariate set. Default is 7.
      }
    
}

\details{
This can be used as an more detailed alternative to labelTopics.
}
\value{
  \item{marginal}{
  A list of matrices, containing the high-probability labels, FREX labels, lift labels,
  and high scoring words.
  }
  \item{K}{
  The number of topics in the STM.
  }
  \item{covnames}{
  Names of the covariate values used in the STM.
  }
  \item{kappa}{Words associated with topics, covariates, and
    topic/covariate interactions.} 
  \item{kappa.m}{Baseline word distribution.}
  \item{n}{
  The n parameter passed by the user to this function; number of words per topic 
  or topic-covariate pair (when covariates are used on the bottom of the model)
  }
  \item{cov.betas}{
  Covariate-specific beta matrices, listing for each covariate a matrix of highest-probability,
  FREX, lift, and high scoring words.  Note that the actual vocabulary has been substituted 
  for word indices.  
  }

}
