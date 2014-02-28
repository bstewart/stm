\name{summary.STM}
\alias{summary.STM}
\alias{print.STM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary Function for the STM objects
}
\description{
Function to report on the contents of STM objects
}
\usage{
\method{summary}{STM}(object, ...)
}
\arguments{
  \item{object}{
An STM object.
}
  \item{\dots}{
Additional arguments affecting the summary
}
}
\details{
Summary prints a short statement about the model and then runs \code{\link{labelTopics}}.
}
