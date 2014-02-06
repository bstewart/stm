\name{plot.estimateEffect}
\alias{plot.estimateEffect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot effect of covariates on topics}
\description{
Plots the marginal effect of a covariate on a set of topics selected by
the user. Before running this, the user should run a function to
simulate the necessary confidence intervals.  See \code{\link{estimateEffect}}.
}
\usage{
\method{plot}{estimateEffect}(x, covariate, model=NULL, topics=x$topics,
                    method="pointestimate",
                    cov.value1=NULL, cov.value2=NULL,
                    npoints=100, nsims=100, ci.level=.95,
                    xlim=NULL, ylim=NULL, ylab="",
                    main="",printlegend=TRUE, labeltype="numbers",
                    n=7,frexw=.5, xlab="",
                    add=FALSE, linecol=NULL, width=25, verbose.labels=TRUE,
                    ...)}
\arguments{
  \item{x}{
    Output of estimateEffect, which calculates simulated betas for plotting
  }
  \item{covariate}{
    String of the name of the main covariate of interest.  All other
    covariates within the formula specified in estimateEffect will be
    kept at their median.
    }
  \item{model}{
    Model output, only necessary if labeltype is "prob", "frex",
    "score", or "lift".
    }
  \item{topics}{
    Topics to plot.}
  \item{method}{
     Method used for plotting.  "pointestimate" estimates mean topic
     proportions for each value of the covariate.  "difference"
     estimates the mean difference in topic proportions for two
     different values of the covariate (cov.value1 and cov.value2
     must be specified).  "continuous" estimates how topic proportions
     vary over the support of a continuous covariate.}
  \item{printlegend}{
    Whether to plot a topic legend in the case of a continuous covariate.
  }
  \item{labeltype}{
    Determines the labeltype for the topics.  The default is "number"
    which prints the topic number.  Other options are "prob", which
    prints the highest probability words, "score", "lift", and "frex",
    from labeltopics (see labeltopics() for more details).  The user can
    also input a vector for labeltype with their own labels.  Labels
    appear in the legend for continous covariates.
  }
  \item{cov.value1}{
    For method "difference", the value or set of values of interest at which to set the covariate.
  }
    \item{cov.value2}{
    For method "difference", the value or set of values which will be
    set as the comparison group.  cov.value1 and cov.value2 must be
    vectors of the same length.
    }
  \item{npoints}{
    Number of unique points to use for simulation along the support of a
    continuous covariate.  For method "continuous" only.
  }
  \item{nsims}{
    Number of simulations for estimation.
    }
  \item{n}{
    Number of words to print if "prob", "score", "lift", or "frex" is chosen. 
  }
  \item{ci.level}{
    Confidence level for confidence intervals.
    }
  \item{frexw}{
    If "frex" labeltype is used, this will be the frex weight.
  }
  \item{add}{
    Logical parameter for whether the line should be added to the plot,
    or a new plot should be drawn.
  }
  \item{linecol}{
    For continuous covariates only.  A vector that specifies the colors
    of the lines within the plot.  If NULL, then colors will be randomly generated.
  }
  \item{verbose.labels}{
    For method "difference" -- verboselabels will specify the comparison
    covariate values
    of the covariate on the plot.
    }
     \item{xlim}{
     Vector of x axis minimum and maximum values.
    }
     \item{ylim}{
     Vector of y axis minimum and maximum values.
    }
         \item{main}{
Character string that is plot title.
}
         \item{xlab}{
Character string that is x axis title.
}
         \item{ylab}{
Character string that is y axis title.
}
         \item{width}{
Character string that is y axis title.
}
    \item{...}{
    Other plotting parameters
  }
}

\details{
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
\examples{
library(lda)
data(poliblog.ratings)
poliblog.ratings <- as.factor(ifelse(poliblog.ratings==-100,
                                     "Liberal",
                                     "Conservative"))
prep <- estimateEffect(1:5 ~ poliblog.ratings, poliblogPrevFit)
plot.estimateEffect(prep, "poliblog.ratings", model=poliblogPrevFit,
method="pointestimate")
plot.estimateEffect(prep, "poliblog.ratings", model=poliblogPrevFit,
method="difference", cov.value1="Liberal", cov.value2="Conservative")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
