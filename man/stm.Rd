\name{stm}
\alias{stm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variational EM for the Structural Topic Model
}
\description{
Estimation of the Structural Topic Model using semi-collapsed variational EM.  The function takes sparse representation of documents, an integer number of topics, and covariates and returns fitted model parameters.  Covariates can be used in the prior for topic \code{prevalence}, in the prior for topical \code{content} or both.
}
\usage{
stm(documents, vocab, K, 
    prevalence, content, data = NULL, 
    init.type = c("LDA", "DMR", "Random"), seed = NULL,
    max.em.its = 100, emtol = 1e-05, 
    verbose = TRUE, ...)
}
\arguments{
   \item{documents}{
  The documents to be modeled.  Object must be a list of with each element corresponding to a document.  Each document is represented as an integer matrix with two rows, and columns equal to the number of unique vocabulary words in the document.  The first row contains the 1-indexed vocabulary entry and the second row contains the number of times that term appears.
 
  This is similar to the format in the \pkg{lda} package except that (following R convention) the vocabulary is indexed from one. Corpora can be imported using the reader function and manipulated using the \code{\link{prepDocuments}}.  Raw texts can be ingested using \code{\link{textProcessor}}.
  }
  \item{vocab}{
  Character vector specifying the words in the corpus in the order of the vocab indices in documents. Each term in the vocabulary index must appear at least once in the documents.  See \code{\link{prepDocuments}} for dropping unused items in the vocabulary.
}
  \item{K}{
  A positive integer (of size 2 or greater) representing the desired number of topics. Additional detail on choosing the number of topics in details.
}
  \item{prevalence}{
  A formula object with no response variable or a matrix containing topic prevalence covariates.  Use \code{s()}, \code{ns()} or \code{bs()} to specify smooth terms. See details for more information.
}
  \item{content}{
  A formula containing a single variable, a factor variable or something which can be coerced to a factor indicating the category of the content variable for each document.
}
  \item{data}{
Dataset which contains prevalence and content covariates.}
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
Convergence tolerance.  EM stops when the relative change in the approximate bound drops below this level.  Defaults to .001\%.
}
  \item{verbose}{
  A logical flag indicating whether information should be printed to the screen.  During the E-step (iteration over documents) a dot will print each time 1\% of the documents are completed.  At the end of each iteration the approximate bound will also be printed.
}
  \item{\dots}{
 Additional options described in details.
}
}
\details{
The main function for estimating a Structural Topic Model (STM).  STM is an admixture with covariates in both mixture components.  Users provide a corpus of documents and a number of topics.  Each word in a document comes from exactly one topic and each document is represented by the proportion of its words that come from each of the K topics.  These proportions are found in the N (number of documents) by K (user specified number of topics) theta matrix.  Each of the K topics are represented as distributions over words.  The K-by-V (number of words in the vocabulary) matrix logbeta contains the natural log of the probability of seeing each word conditional on the topic.  

The most important user input in parametric topic models is the number of topics.  There is no right answer to the appropriate number of topics.  More topics will give more fine-grained representations of the data at the potential cost of being less precisely estimated.  The number must be at least 2 which is equivalent to unidimensional scaling model.  For short corpora focused on very specific subject matter (such as survey experiments) 3-5 topics is a useful starting range.  For small corpora (a few hundred to a few thousand) 5-20 topics is a good place to start.  Beyond these rough guidelines it is application specific.  Previous applications in political science with medium sized corpora (10k to 100k documents) have found 50-60 topics to work well.  For larger corpora 100 topics is a useful default size.  Of course, your mileage may vary.

The model for topical prevalence includes covariates which the analyst believes may influence the frequency with which a topic is discussed.  This is specified as a formula which can contain smooth terms using splines or by using the function \code{\link{s}}.  The response portion of the formula should be left blank.  See the examples.

The topical convent covariates are those which affect the way in which a topic is discussed. As currently implemented this must be a single variable which defines a discrete partition of the dataset (each document is in one and only one group).  We may relax this in the future.  While including more covariates in topical prevalence will rarely affect the speed of the model, including additional levels of the content covariates can make the model much slower to converge.  This is due to the model operating in the much higher dimensional space of words in dictionary (which tend to be in the thousands) as opposed to topics.  

Additional options can be specified as follows:
\enumerate{
  \item \code{LDAbeta} is a logical that defaults to \code{TRUE} when there are no content covariates.  When set to \code{FALSE} the model performs SAGE style topic updates (sparse deviations from a baseline).
  \item \code{interactions} is a logical that defaults to \code{TRUE}.  This automatically includes interactions between content covariates and the latent topics.  Setting it to \code{FALSE} reduces to a model with no interactive effects.
  \item \code{gammamode} sets the prior estimation method for the prevalence covariate model.  The default \code{Pooled} options uses Normal prior distributions with a topic-level pooled variance which is given a broad gamma hyperprior.  The alternative \code{GL} uses Taddy's Gamma-Lasso which implies Laplace prior distributions with element-wise variances given gamma priors.  The \code{GL} option is useful when sparsity in the covariate effects is desired.
  \item \code{sigmaprior} a scalar between 0 and 1 which defaults to 0.  This sets the strength of regularization towards a diagonalized covariance matrix.  Setting the value above 0 can be useful if topics are becoming too highly correlated
  \item \code{taumode} sets the prior estimation for the content covariate coefficients.  The default is \code{Jeffreys} and uses the improper Jeffreys prior.  The option \code{GL} uses the Gamma Lasso.  Both options are sparsity promoting.
  \item \code{reportevery} sets when top words are printed to the screen.  Defaults to 5.
  \item \code{keepHistory} a logical defaulting to \code{FALSE} which indicates whether the history of parameter values at each level should be saved.  Note that this can be extremely memory intensive.
}
}
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
  Evaluation. http://scholar.harvard.edu/files/bstewart/files/stmnips2013.pdf


  Roberts, M., Stewart, B., Tingley, D., Lucas, C., Leder-Luis, J.,
  Gadarian, S., Albertson, B., Albertson, B. and Rand,
  D. (Forthcoming). "Structural topic models for open ended survey
  responses." American Journal of Political Science http://scholar.harvard.edu/files/dtingley/files/topicmodelsopenendedexperiments.pdf
}

\seealso{
\code{\link{prepDocuments}}
\code{\link{labelTopics}}
\code{\link{estimateEffect}}
}
\examples{
library(lda)

#Let's Start By Loading a Blog Corpus in the LDA Package
data(poliblog.documents) # these are the word counts
data(poliblog.vocab)     # these are the words in the vocabulary
data(poliblog.ratings)

#The rating is in a rather unintuitive format so let's recode it to a factor variable
poliblog.ratings <- as.factor(ifelse(poliblog.ratings==-100, "Liberal", "Conservative"))

#Let's see what we have
length(poliblog.documents) #there are 773 documents
length(poliblog.vocab)     #and 1290 items in the vocab

#The lda package uses a slightly different format than us.  
# specifically they index words from 0 using the convention in C.
# So let's convert it over to a more R-like 1-index.
out <- prepDocuments(poliblog.documents, poliblog.vocab)
#the prepDocuments() function actually does a lot of other things but for now we are just
# using it to deal with the zero-indexing

#now we reassign the documents object
poliblog.documents <- out$documents

# stm() is the main function for estimating models
# the simplest possible specification requires documents, vocab and K the number of topics.
# without any covariates the model is closely related to the standard LDA model
# (its actually a correlated topic model see Blei and Lafferty 2006)
mod.out <- stm(poliblog.documents,poliblog.vocab,K=5, max.em.its=5)

#we can quickly look at the results 
summary(mod.out)
#this shows us the top words ranked according to four different criteria.

#We can also include our ratings variable as a covariate in the topic prevalence
# portion of the model.  
mod.out <- stm(poliblog.documents,poliblog.vocab,K=5, prevalence=~poliblog.ratings, max.em.its=5)
summary(mod.out)

#We can also include a covariate in the content equation.  It need not be the same 
\dontrun{
mod.out2 <- stm(poliblog.documents,poliblog.vocab,K=5,  prevalence=~poliblog.ratings,
                content=~poliblog.ratings, max.em.its=1)
summary(mod.out2)
}

}