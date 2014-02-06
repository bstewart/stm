\name{topicQuality}
\alias{topicQuality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots semantic coherence and exclusivity for each topic.}
\description{
  Plots semantic coherence and exclusivity for each topic. 
  Does not support models.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
topicQuality(model, documents, xlab="Semantic Coherence",ylab="Exclusivity", labels=1:ncol(model$theta), M=10,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    Output from stm, or a selected model from selectModel.
  }
  \item{documents}{
    Vector containing documents used.
  }
    \item{labels}{
Vector of number corresponding to topic numbers.
}
  \item{M}{
    Number of words to use in semantic coherence and exclusivity calculations
  }
           \item{xlab}{
Character string that is x axis title. This should be semantic coherence.
}
         \item{ylab}{
Character string that is y axis title. This should be exclusivity.
}
  \item{...}{
    Other plotting parameters from igraph.
  }
}

\details{
  Each model has semantic coherence and exclusivity values associated
  with each topic.  This function plots these values and labels each with its topic number.

%%  ~~ If necessary, more details than the description above ~~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
