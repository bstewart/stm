\name{stm-package}
\alias{stm-package}
\docType{package}
\title{
Structural Topic Model
}
\description{
This package implements the Structural Topic Model and eventually other
meta-data models.  The packages is still in beta so small details of the
function syntax may change.  \code{\link{stm}} is the main function of the package that estimates the
Structural Topic Model.  Data import can be done via the read functions
or if from raw texts the textProcessor function.  \code{\link{prepDocuments}} covers
some basic format conversions and will do other sorts of data
manipulation.  \code{\link{labelTopics}} will show the topics fitted by the model.
\code{\link{estimateEffect}}, \code{\link{plot.estimateEffect}}
will show relationships between
topics and covariates within the model.   Various other plotting
functions are included that allow the user to plot prevalence of topics,
correlations between topics, and other types of summaries of the model.
}
\details{
\tabular{ll}{
Package: \tab stm\cr
Type: \tab Package\cr
Version: \tab 0.05.03\cr
Date: \tab 2014-02-13\cr
License: \tab MIT \cr
}
}
\author{
Author: Margaret E. Roberts, Brandon M. Stewart and Dustin Tingley

Maintainer: Brandon Stewart <bstewart@fas.harvard.edu>
}
\references{
Roberts, M., Stewart, B., Tingley, D., and Airoldi, E. (2013) "The
structural topic model and applied social science." In Advances in
Neural Information Processing Systems Workshop on Topic Models:
Computation, Application, and Evaluation. http://scholar.harvard.edu/files/bstewart/files/stmnips2013.pdf

Roberts, M., Stewart, B., Tingley, D., Lucas, C., Leder-Luis, J.,
Gadarian, S., Albertson, B., Albertson, B. and Rand,
D. (Forthcoming). "Structural topic models for open ended survey
responses." American Journal of Political Science http://scholar.harvard.edu/files/dtingley/files/topicmodelsopenendedexperiments.pdf

}
\keyword{ package }
\seealso{
\code{\link{stm}}
}

