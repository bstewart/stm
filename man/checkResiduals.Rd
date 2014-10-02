\name{checkResiduals}
\alias{checkResiduals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Residual dispersion test for topic number
  }
\description{
Computes the multinomial dispersion of the STM residuals as in Taddy (2012)
}
\usage{
checkResiduals(stmobj, documents, tol=.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stmobj}{
     an \code{STM} model object for which to compute residuals.
      }
  \item{documents}{
	 the documents corresponding to \code{stmobj} as in \code{\link{stm}}.
      }
  \item{tol}{
  the tolerance parameter for calculating the degrees of freedom.  Defaults to 1/100 as in Taddy(2012)
  }}
    
\details{
This function implements the residual-based diagnostic method of Taddy (2012).  The basic idea is that when the model is correctly specified the multinomial likelihood implies a dispersion of the residuals: \eqn{\sigma^2=1}.  If we calculate the sample dispersion and the value is greater than one, this implies that the number of topics is set too low, because the latent topics are not able to account for the overdispersion. In practice this can be a very demanding criterion, especially if the documents are long.  However, when coupled with other tools it can provide a valuable perspective on model fit.

Further details are available in the referenced paper, but broadly speaking the dispersion is derived from the mean of the squared adjusted residuals.  We get the sample dispersion by dividing by the degrees of freedom parameter.  In estimating the degrees of freedom, we follow Taddy (2012) in approximating the parameter \eqn{\hat{N}} by the number of expected counts exceeding a tolerance parameter.  The default value of 1/100 given in the Taddy paper can be changed by setting the \code{tol} argument.  

The function returns the estimated sample dispersion (which equals 1 under the data generating process) and the p-value of a chi-squared test where the null hypothesis is that \eqn{\sigma^2=1} vs the alternative \eqn{\sigma^2 >1}. As Taddy notes and we echo, rejection of the null 'provides a very rough measure for evidence in favor of a larger number of topics.'
}
\references{
     Taddy, M. 'On Estimation and Selection for Topic Models'. AISTATS 2012, JMLR W&CP 22
}
\examples{
\dontrun{
#An example using the Gadarian data.  From Raw text to fitted model.
temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
meta<-temp$meta
vocab<-temp$vocab
docs<-temp$documents
out <- prepDocuments(docs, vocab, meta)
docs<-out$documents
vocab<-out$vocab
meta <-out$meta
set.seed(02138)
mod.out <- stm(docs, vocab, 3, prevalence=~treatment + s(pid_rep), data=meta)
checkResiduals(mod.out, docs)
}
}
