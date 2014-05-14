\name{permutationTest}
\alias{permutationTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutation test of a binary covariate.
}
\description{
Run a permutation test where a binary treatment variable is randomly permuted and topic model is reestimated.
}
\usage{
permutationTest(formula, stmobj, treatment, 
                nruns=100, 
                documents, vocab, data, seed=NULL,
                stmverbose=TRUE,uncertainty="Global")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the prevalence component of the \code{stm} model and the \code{estimateEffect} call.  This formula must contain at least one binary covariate (specified through the argument \code{treatment}) but it can contain other terms as well.  If the binary covariate is interacted with additional variables the estimated quantity of interest is the effect when those additional variables are set to 0.
}
  \item{stmobj}{
    Model output from a single run of \code{stm} which contains the reference effect.
}
  \item{treatment}{
  A character string containing treatment id as used in the formula of the stmobj.  This is the variable which is randomly permuted.
}
  \item{nruns}{
  Number of total models to fit (including the original model).
  }
  \item{documents}{
  The documents used in the stmobj model.
  }
  \item{vocab}{
  The vocab used in the stmobj model.
  }
  \item{data}{
  The data used in the stmobj model.
  }
  \item{seed}{
  Optionally a seed with which to replicate the result.  As in \code{\link{stm}} the seed is automatically saved and returned as part of the object.  Passing the seed here will replicate the previous run.
  }
  \item{stmverbose}{
  Should the stm model be run with \code{verbose=TRUE}.  Turning this to \code{FALSE} will surpress only the model specific printing.  An update on which model is being run will still print to the screen.
  }
  \item{uncertainty}{
Which procedure should be used to approximate the measurement uncertainty in the topic proportions.  See details for more information.  Defaults to the Global approximation.
}
}
\details{
This function takes a single binary covariate and runs a permutation test where, rather than using the true assignment, the covariate is randomly drawn with probability equal to its empirical probability in the data.
After each shuffle of the covariate the same STM model is estimated at different starting values using the same initialization procedure as the original model, and the effect of the covariate across topics is calculated.

Next the function records two quantities of interest across this set of "runs" of the model. The first records the absolute maximum effect of the permuted covariate across all topics.

The second records the effect of the (permuted) covariate on the topic in each additional stm run which is estimated to be the topic closest to the topic of interest (specified in \code{\link{plot.STMpermute}}) from the original stm model. Uncertainty can be calculated using the standard options in \code{\link{estimateEffect}}.
}
\value{
\item{ref}{A list of K elements containing the quantiles of the estimated effect for the reference model.}
\item{permute}{A list where each element is an aligned model parameter summary}
\item{variable}{The variable id that was permuted.}
\item{seed}{The seed for the stm model.}
}
\seealso{
\code{\link{plot.STMpermute}}
}
\examples{
\dontrun{
temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
documents <- out$documents
vocab <- out$vocab
meta <- out$meta
set.seed(02138)
mod.out <- stm(documents, vocab, 3, prevalence=~treatment + s(pid_rep), data=meta)
summary(mod.out)
prep <- estimateEffect(1:3 ~ treatment + s(pid_rep), mod.out, meta)
plot.estimateEffect(prep, "treatment", model=mod.out,
                    method="difference",cov.value1=1,cov.value2=0)
test <- permutationTest(formula=~ treatment + s(pid_rep), stmobj=mod.out, 
                        treatment="treatment", nruns=25, documents=documents,
                        vocab=vocab,data=meta, stmverbose=FALSE)
plot.STMpermute(test,2, xlab="Effect", ylab="Model Index", main="Topic 2 Placebo Test")
}
}
