\name{factorCheck}
\alias{factorCheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Control for topics loading exclusively on a single covariate level
}
\description{
Controls for the presence of aberrations in a STM model with a discrete covariate. Computes topical proportions by covariate level, and outputs a warning if any value is above a user-specified threshold, signifying an abnormal concentration of a topic's mass on a single covariate level.
}
\usage{
factorCheck(stmobject, covariate=NULL, tolerance=0.01, reporting=0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stmobject}{
The STM model object.
}

\item{covariate}{
A vector containing the relevant covariate levels.
}

\item{tolerance}{
The minimum topical proportion (by document) to consider a document as containing a particular topic.
}

\item{reporting}{
The minimum threshold for topical proportion (by covariate level) to report an abberration.
}
}
\details{
This function computes topical proportion by covariate level using the following algorithm. First, the theta matrix in the STM object (topical proportions by document) is subsetted by covariate levels. Second, indicators are computed specifying whether the entries in the subsetted matrix are higher than the tolerance value. Third, the indicators are aggregated by topic, summing over the columns of the subsetted matrix. Fourth, the aggregate value is divided by the total number of observations above the tolerance value in the corresponding column of the non-subsetted theta matrix, to obtain the topical proportion by covariate level.
}
\value{
The function returns a list of two items:
\item{topical.proportions}{
 A matrix containing topical proportions by covariate levels. Rows correspond to levels, columns to topics.
}

\item{total.errors}{
The number of matrix entries above the given reporting threshold.
}
}

\author{
Antonio Coppola, Harvard University: <acoppola@college.harvard.edu>
}


\examples{
## Perform test on Gadarian and Albertson data
diagnostic <- factorCheck(gadarianFit, gadarian$treatment)
}

