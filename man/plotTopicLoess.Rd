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
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.estimateEffect}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
