\name{stm}
\alias{stm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variational EM for the Structural Topic Model
}
\description{
Estimation of the Structural Topic Model using semi-collapsed variational EM.  The function takes sparse representation of documents, an integer number of topics, and covariates and returns fitted model parameters.  Covariates can be used in the prior for topic \code{prevalence}, in the prior for topical \code{content} or both.  See an overview of functions in the package here: \code{\link{stm-package}}
}
\usage{
stm(documents, vocab, K, 
    prevalence, content, data=NULL,
    init.type=c("LDA", "DMR","Random"), seed=NULL, 
    max.em.its=100, emtol=1e-5,
    verbose=TRUE, reportevery=5, keepHistory=FALSE,  
    LDAbeta=TRUE, interactions=TRUE,
    gamma.prior=c("Pooled", "L1"), sigma.prior=0,
    kappa.prior=c("Jeffreys", "L1"), control=list()) 
}
\arguments{
   \item{documents}{
  The documents to be modeled.  Object must be a list of with each element corresponding to a document.  Each document is represented as an integer matrix with two rows, and columns equal to the number of unique vocabulary words in the document.  The first row contains the 1-indexed vocabulary entry and the second row contains the number of times that term appears.
 
  This is similar to the format in the \code{\link[lda]{lda}} package except that (following R convention) the vocabulary is indexed from one. Corpora can be imported using the reader function and manipulated using the \code{\link{prepDocuments}}.  Raw texts can be ingested using \code{\link{textProcessor}}.
  }
  \item{vocab}{
  Character vector specifying the words in the corpus in the order of the vocab indices in documents. Each term in the vocabulary index must appear at least once in the documents.  See \code{\link{prepDocuments}} for dropping unused items in the vocabulary.
}
  \item{K}{
  A positive integer (of size 2 or greater) representing the desired number of topics. Additional detail on choosing the number of topics in details.
}
  \item{prevalence}{
  A formula object with no response variable or a matrix containing topic prevalence covariates.  Use \code{\link{s}}, \code{\link[splines]{ns}} or \code{\link[splines]{bs}} to specify smooth terms. See details for more information.
}
  \item{content}{
  A formula containing a single variable, a factor variable or something which can be coerced to a factor indicating the category of the content variable for each document.
}
  \item{data}{
an optional data frame containing the prevalence and/or content covariates.  If unspecified the variables are taken from the active environment.}
  \item{init.type}{
The method of initialization.  Must be either Latent Dirichlet Allocation (LDA), Dirichlet Multinomial Regression Topic Model (DMR), or a random initialization.  If you want to replicate a previous result, see the argument \code{seed}.  
}
  \item{seed}{
Seed for the random number generator. \code{stm} saves the seed it uses on every run so that any result can be exactly reproduced.  When attempting to reproduce a result with that seed, it should be specified here.
}
  \item{max.em.its}{
  The maximum number of EM iterations.  If convergence has not been met at this point, a message will be printed.
}
  \item{emtol}{
Convergence tolerance.  EM stops when the relative change in the approximate bound drops below this level.  Defaults to .001.
}
  \item{verbose}{
  A logical flag indicating whether information should be printed to the screen.  During the E-step (iteration over documents) a dot will print each time 1\% of the documents are completed.  At the end of each iteration the approximate bound will also be printed.
}
  \item{reportevery}{
  An integer determining the intervals at which labels are printed to the screen during fitting.  Defaults to every 5 iterations.
}
 \item{keepHistory}{
 Logical indicating whether the history should be saved at each iteration.  Defaults to \code{FALSE}.  Note that the model parameters are extremely memory intensive so use with care.
 }
 \item{LDAbeta}{
 a logical that defaults to \code{TRUE} when there are no content covariates.  When set to \code{FALSE} the model performs SAGE style topic updates (sparse deviations from a baseline).
 }
 \item{interactions}{
 a logical that defaults to \code{TRUE}.  This automatically includes interactions between content covariates and the latent topics.  Setting it to \code{FALSE} reduces to a model with no interactive effects.
 }
 \item{gamma.prior}{
 sets the prior estimation method for the prevalence covariate model.  The default \code{Pooled} options uses Normal prior distributions with a topic-level pooled variance which is given a broad gamma hyperprior.  The alternative \code{L1} uses \code{glmnet} to estimate a grouped penalty between L1-L2.  See details below.
 }
 \item{sigma.prior}{
 a scalar between 0 and 1 which defaults to 0.  This sets the strength of regularization towards a diagonalized covariance matrix.  Setting the value above 0 can be useful if topics are becoming too highly correlated
 }
 \item{kappa.prior}{
 sets the prior estimation for the content covariate coefficients.  The default is \code{Jeffreys} and uses a scale mixture of Normals with an improper Jeffreys prior.  The option \code{L1} uses \code{glmnet} to estimate with a penalty between L1 and L2,  See details for more.
 }
 \item{control}{
 a list of additional parameters control portions of the optimization.  See details.
 }
}
\details{
The main function for estimating a Structural Topic Model (STM).  STM is an admixture with covariates in both mixture components.  Users provide a corpus of documents and a number of topics.  Each word in a document comes from exactly one topic and each document is represented by the proportion of its words that come from each of the K topics.  These proportions are found in the N (number of documents) by K (user specified number of topics) theta matrix.  Each of the K topics are represented as distributions over words.  The K-by-V (number of words in the vocabulary) matrix logbeta contains the natural log of the probability of seeing each word conditional on the topic.  

The most important user input in parametric topic models is the number of topics.  There is no right answer to the appropriate number of topics.  More topics will give more fine-grained representations of the data at the potential cost of being less precisely estimated.  The number must be at least 2 which is equivalent to unidimensional scaling model.  For short corpora focused on very specific subject matter (such as survey experiments) 3-5 topics is a useful starting range.  For small corpora (a few hundred to a few thousand) 5-20 topics is a good place to start.  Beyond these rough guidelines it is application specific.  Previous applications in political science with medium sized corpora (10k to 100k documents) have found 50-60 topics to work well.  For larger corpora 100 topics is a useful default size.  Of course, your mileage may vary.

The model for topical prevalence includes covariates which the analyst believes may influence the frequency with which a topic is discussed.  This is specified as a formula which can contain smooth terms using splines or by using the function \code{\link{s}}.  The response portion of the formula should be left blank.  See the examples.

The topical convent covariates are those which affect the way in which a topic is discussed. As currently implemented this must be a single variable which defines a discrete partition of the dataset (each document is in one and only one group).  We may relax this in the future.  While including more covariates in topical prevalence will rarely affect the speed of the model, including additional levels of the content covariates can make the model much slower to converge.  This is due to the model operating in the much higher dimensional space of words in dictionary (which tend to be in the thousands) as opposed to topics.  

In addition to the default priors for prevalence and content, we also make use of the \code{glmnet} package to allow for penalties between the L1 and L2 norm.  In these settings we estimate a regularization path and then select the optimal shrinkage parameter using a user-tuneable information criterion.  By default selecting the \code{L1} option will apply the L1 penalty selecting the optimal shrinkage parameter using AIC. The defaults have been specifically tuned for the STM but almost all the relevant arguments can be changed through the control structure below.  Changing the \code{gamma.enet} and \code{kappa.enet} parameters allow the user to choose a mix between the L1 and L2 norms.  When set to 1 (as by default) this is the lasso penalty, when set to 0 its the ridge penalty.  Any value in between is a mixture called the elastic net.  

The \code{control} argument is a list with named components which can be used to specify numerous additional computational details.  Valid components include:
  \describe{
    \item{\code{tau.maxit}}{Controls the maximum number of iterations when estimating the prior for content covariates.  When the mode is \code{Jeffreys} as by default, estimation proceeds by iterating between the kappa vector corresponding to a particular topic and the associated variance tau before moving on to the next parameter vector. this controls the maximum number of iterations. It defaults to \code{NULL} effectively enforcing convergence.  When the mode is \code{L1} this sets the maximum number of passes in the coordinate descent algorithm and defaults to 1e8.}
    \item{\code{tau.tol}}{Sets the convergence tolerance in the optimization for content covariates.  When the mode is \code{Jeffreys} this sets the convergence tolerance in the iteration between the kappa vector and variances tau and defaults to 1e-5.  With \code{L1} it defaults to 1e-6.}
    \item{\code{kappa.mstepmaxit}}{When the mode for content covariate estimation is \code{Jeffreys} this controls the maximum number of passes through the sequence of kappa vectors.  It defaults to 3.  It has no role under \code{L1}- see \code{tau.maxit} option instead.}
    \item{\code{kappa.msteptol}}{WHen the mode for content covariate estimation is \code{Jeffreys} this controls the tolerance for convergence (measured by the L1 norm) for the entire M-step.  It is set to .99 by default.  This has no role under mode \code{L1}- see \code{tau.tol} option instead.}
    \item{\code{wordconverge.num}}{provides an alternative convergence metric based on the number of iterations without change in the top n most probable words per topic.  This sets the number n.  This is 20 by default (although it is deactivated by default as well- see the next argument). 
    }
    \item{\code{fixedintercept}}{a logical indicating whether in content covariate models the intercept should be fixed to the background distribution.  TRUE by default. This only applies when kappa.prior is set to L1.  If FALSE the intercept is estimated from the data without penalty.  In practice estimated intercepts often push term probabilities to zero, resulting in topics that look more like those in a Dirichlet model- that is, most terms have approximately zero probability with some terms with high probability.}
    \item{\code{kappa.enet}}{When using the L1 mode for content covariates this controls the elastic net mixing parameter.  See the argument \code{alpha} in \code{glmnet}.  Value must be between 1 and 0 where 1 is the lasso penalty (the default) and 0 is the ridge penalty.  The closer the parameter is to zero the less sparse the solution will tend to be.}
    \item{\code{gamma.enet}}{Controls the elastic net mixing parameter for the prevalence covariates.  See above for a description.}
    \item{\code{nlambda}}{Controls the length of the regularization path when using L1 mode for content covariates.  Defaults to 500.  Note that glmnet relies heavily on warm starts and so a high number will often (counter-intuitively) be less costly than a low number.  We have chosen a higher default here than the default in the glmnet package and we don't recommend changing it.}
    \item{\code{lambda.min.ratio}}{For L1 mode content covariates this controls the explored path of regularization values.  This defaults to .0001.  Setting higher numbers will result in more sparse solutions.  This is here primarily for dealing with convergence issues, if you want to favor selection of sparser solutions see the next argument.}
    \item{\code{ic.k}}{For L1 mode content covariates this controls the selection of the regularization parameter.  We use a generic information criterion which penalizes complexity by the parameter \code{ic.k}.  When set to 2 (as by default) this results in AIC.  When set to log(n) (where n is the total number of words in the corpus) this is equivalent to BIC.  Larger numbers will express a preference for sparser (simpler) models.}
}}
\value{
An object of class STM
\item{mu }{The corpus mean of topic prevalence and coefficients}
\item{sigma }{Covariance matrix}
\item{beta }{List containing the log of the word probabilities for each topic.}
\item{settings }{The settings file.  The Seed object will always contain the seed which can be fed as an argument to recover the model.}
\item{vocab }{The vocabulary vector used.}
\item{convergence }{list of convergence elements including the value of the approximate bound on the marginal likelihood at each step.}
\item{theta}{ Number of Documents by Number of Topics matrix of topic proportions.}
\item{eta}{Matrix of means for the variational distribution of the multivariate normal latent variables used to calculate theta.}
\item{history }{If keepHistory=TRUE the history of model parameters at each step.}
}

\references{
  Roberts, M., Stewart, B., Tingley, D., and Airoldi, E. (2013) "The
  structural topic model and applied social science." In Advances in
  Neural Information Processing Systems Workshop on Topic Models:
  Computation, Application, and
  Evaluation. http://goo.gl/uHkXAQ

  Roberts M., Stewart, B. and Airoldi, E. (2014) "Structural Topic Models" 
  
  Roberts, M., Stewart, B., Tingley, D., Lucas, C., Leder-Luis, J.,
  Gadarian, S., Albertson, B., Albertson, B. and Rand,
  D. (Forthcoming). "Structural topic models for open ended survey
  responses." American Journal of Political Science http://goo.gl/0x0tHJ
}

\seealso{
\code{\link{prepDocuments}}
\code{\link{labelTopics}}
\code{\link{estimateEffect}}
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


}
}
