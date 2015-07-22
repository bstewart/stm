\name{topicLasso}
\alias{topicLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot predictions using topics
}
\description{
Use the \pkg{glmnet} package to plot LASSO based estimates of relationship between an arbitrary dependent variable with topics and additional variables as predictors.  This function is experimental (see below).
}
\usage{
topicLasso(formula, data, stmobj=NULL, subset=NULL,
                      omit.var=NULL,
                      family="gaussian",  main="Topic Effects on Outcome",
                      xlab=expression("Lower Outcome                  Higher Outcome"),
                      labeltype=c("prob", "frex", "lift", "score"),seed=02138,
                      xlim=c(-4,4), standardize=FALSE, nfolds=20, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Formula specifying the dependent variable and additional variables to included in the LASSO beyond the topics present in the stmobj.  Just pass a 1 on the right-hand side in order to run without additional controls.
}
\item{data}{
Data file containing the dependent variable. Typically will be the metadata file used in the stm analysis. It must have a number of rows equal to the number of documents in the stmobj.
}
\item{stmobj}{
The STM object, and output from the \code{stm} function.
}
\item{subset}{
A logical statement that will be used to subset the corpus.
}
\item{omit.var}{
Pass a character vector of variable names to be excluded from the plot.  Note this does not exclude them from the calculation, only the plot.
}
\item{family}{
The family parameter used in \code{\link{glmnet}}.  See explanation there.  Defaults to "gaussian"
}
\item{main}{
Character string for the main title.
}
\item{xlab}{
Character string giving an x-axis label.
}
\item{labeltype}{
Type of example words to use in labelling each topic. See \code{\link{labelTopics}}. Defaults to "prob".
}
\item{seed}{
The random seed for replication of the cross-validation samples.
}
\item{xlim}{
Width of the x-axis.
}
\item{standardize}{
Whether to standardize variables. Default is FALSE, which is different from the \pkg{glmnet} default because the topics are already standarized.  Note that glmnet standardizes the variables by default but then projects them back to their original scales before reporting coefficients.
}
\item{nfolds}{
the number of cross-validation folds.  Defaults to 20.
}
\item{...}{
Additional arguments to be passed to glmnet.  This can be useful for addressing convergence problems.
}
}
\details{
This function is used for estimating the most important topical predictors of an arbitrary outcome.  The idea is to run an L1 regularized regression using \code{\link{cv.glmnet}} in the \pkg{glmnet} package where the document-level dependent variable is chosen by the user and the predictors are the document-topic proportions in the \code{\link{stm}} model along with any other variables of interest.

The function uses cross-validation to choose the regularization parameter and generates a plot of which loadings were the most influential in predicting the outcome.  It also invisibly returns the glmnet model so that it can be used for prediction.

NOTE: This function is still very experimental and may have stability issues.  If stability issues are encountered see the documentation in \pkg{glmnet} for arguments that can be passed to improve convergence.  Also, it is unlikely to work well with multivariate gaussian or multinomial families.
}
\references{
Friedman, Jerome, Trevor Hastie, and Rob Tibshirani. "Regularization paths for generalized linear models via coordinate descent." Journal of statistical software 33.1 (2010): 1.
}
\seealso{
\pkg{glmnet}
}

\examples{
\dontrun{
#Load the poliblog data
data(poliblog5k)
#estimate a model with 50 topics
stm1 <- stm(poliblog5k.docs, poliblog5k.voc, 50,
            prevalence=~rating + blog, data=poliblog5k.meta,
            init.type="Spectral")
#make a plot of the topics most predictive of "rating"
out <- topicLasso(rating ~ 1, family="binomial", data=poliblog5k.meta,stmobj=stm1)
#generate some in-sample predictions
pred <- predict(out, newx=stm1$theta,type="class")
#check the accuracy of the predictions
table(pred, poliblog5k.meta$rating)
}
}
