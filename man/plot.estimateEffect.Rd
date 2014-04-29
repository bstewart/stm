\name{plot.estimateEffect}
\alias{plot.estimateEffect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot effect of covariates on topics}
\description{
Plots the effect of a covariate on a set of topics selected by
the user. Different effect types available depending on type of covariate. Before running this, the user should run a function to
simulate the necessary confidence intervals.  See \code{\link{estimateEffect}}.
}
\usage{
\method{plot}{estimateEffect}(x, covariate, model=NULL, topics=x$topics,
                          method="pointestimate",
                          cov.value1=NULL, cov.value2=NULL,int.value=NULL, 
                          npoints=100, nsims=100, ci.level=.95, 
                          xlim=NULL, ylim=NULL, ylab=NULL, 
                          main="",printlegend=TRUE, 
                          labeltype="numbers",n=7,frexw=.5, 
                          xlab="", add=FALSE, linecol=NULL, 
                          width=25, verbose.labels=TRUE, 
                          family=NULL, custom.labels=NULL,...)}
\arguments{
  \item{x}{
    Output of estimateEffect, which calculates simulated betas for plotting
  }
  \item{covariate}{
    String of the name of the main covariate of interest. Must be enclosed in quotes.  All other
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
       \item{cov.value1}{
    For method "difference", the value or set of values of interest at which to set the covariate. In the case of calculating a treatment/control contrast, set the treatment to cov.value1.
  }
    \item{cov.value2}{
    For method "difference", the value or set of values which will be
    set as the comparison group.  cov.value1 and cov.value2 must be
    vectors of the same length.
  }
  \item{int.value}{
    For methods "pointestimate" and "difference", when two binary terms
    are interacted and one variable in the interaction is the covariate
    of interest, the user can specify the value of the interaction term.
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
  \item{printlegend}{
    Whether to plot a topic legend in the case of a continuous covariate.
  }
  \item{labeltype}{
    Determines the labeltype for the topics.  The default is "number"
    which prints the topic number.  Other options are "prob", which
    prints the highest probability words, "score", "lift", and "frex",
    from labeltopics (see labeltopics() for more details).  The user can
    also select "custom" for custom labels, which should be inputted
    under custom.labels. Labels
    appear in the legend for continous covariates.
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
\item{custom.labels}{
  A vector of custom.labels if labeltype is equal to "custom".
  }
    \item{...}{
    Other plotting parameters
  }
}
\examples{
\dontrun{
prep <- estimateEffect(1:3 ~ treatment, gadarianFit, gadarian)
plot.estimateEffect(prep, "treatment", model=gadarianFit,
method="pointestimate")
plot.estimateEffect(prep, "treatment", model=gadarianFit,
method="pointestimate")
plot.estimateEffect(prep, "treatment", model=gadarianFit,
method="difference",cov.value1=1,cov.value2=0)

#If the covariate were a binary factor, the factor labels can be used to
 specify the values of cov.value1 (e.g., cov.value1="treat"). String
 variables must be turned to factors prior to plotting.

#Example of binary times binary interaction
gadarian$binaryvar <- sample(c(0,1), nrow(gadarian), replace=T)
temp <- textProcessor(gadarian$open.ended.response,metadata=gadarian)
out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
stm1 <- stm(out$documents, out$vocab, 3, prevalence=~treatment*binaryvar,
 data=gadarian)
prep <- estimateEffect(c(2) ~ treatment*binaryvar, stmobj=stm1,
metadata=gadarian)

par(mfrow=c(1,2))
plot.estimateEffect(prep, "binaryvar", method="pointestimate",
cov.value1=1, cov.value2=0, xlim=c(-1,1), int.value=1)
plot.estimateEffect(prep, "binaryvar", method="pointestimate",
cov.value1=1, cov.value2=0, xlim=c(-1,1), int.value=0)
}
}
