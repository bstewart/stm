\name{plotTopicLoess}
\alias{plotTopicLoess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot some effects with loess}
\description{
  Plots a loess line of the topic proportions on a covariate inputted by
  the user. This allows for a more flexible functional form for the relationship. 
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotTopicLoess(model, topics, covariate, span=1.5, level=.95,
               main="", xlab="Covariate", ylab="Topic Proportions")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    An STM model object
    }
  \item{topics}{
    Vector of topic numbers to plot by the covariate. E.g., c(1,2,3) would plot lines for topics 1,2,3.
  }
  \item{covariate}{
    Covariate vector by which to plot topic proportions.
  }
  \item{span}{
    loess span parameter.  See \code{\link{loess}}
  }
  \item{level}{
    Desired coverage for confidence intervals
  }
  \item{main}{
    Title of the plot, default is ""
  }
  \item{xlab}{
    X-label, default is "Covariate"
  }
  \item{ylab}{
    Y-label, default is "Topic Proportions"
    }
}

\details{
This function is considerably less developed than \code{\link{plot.estimateEffect}} and we recommend using that function with splines and high degrees of freedom where possible.  Computes standard errors through the method of composition as in \code{\link{estimateEffect}}.
}
\seealso{
\code{\link{plot.estimateEffect}}}

