\name{s}
\alias{s}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a B-spline Basis Function
}
\description{
This is a simple wrapper around the \code{\link[splines]{bs}} function in the splines package.  It will default to a spline with 10 degrees of freedom.}
\usage{
s(x, df, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the predictor value.
}
  \item{df}{
degrees of freedom.  Defaults to the minimum of 10 or one minus the number of unique values in x.
}
  \item{\dots}{
Arguments passed to the \code{\link[splines]{bs}} function.
}
}
\details{
This is a simple wrapper written as users may find it easier to simply type \code{\link{s}} rather than selecting parameters for a spline.
}
\value{
a predictor matrix of the basis functions.
}
\seealso{
\code{\link[splines]{bs}}
\code{\link[splines]{ns}}
}

